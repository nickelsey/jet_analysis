// basic_dijet_imbalance.cxx

#include <iostream>
#include <string>
#include <fstream>

#include "util/arg_helper.hh"
#include "util/trigger_lookup.hh"

#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TH2.h"
#include "TH1.h"

#include "TStarJetPicoReader.h"
#include "TStarJetPicoEvent.h"
#include "TStarJetPicoEventHeader.h"
#include "TStarJetPicoEventCuts.h"
#include "TStarJetPicoPrimaryTrack.h"
#include "TStarJetPicoTower.h"
#include "TStarJetPicoTrackCuts.h"
#include "TStarJetPicoTowerCuts.h"
#include "TStarJetVectorContainer.h"
#include "TStarJetVector.h"
#include "TStarJetPicoTriggerInfo.h"
#include "TStarJetPicoUtils.h"

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/ClusterSequencePassiveArea.hh"
#include "fastjet/ClusterSequenceActiveArea.hh"
#include "fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/FunctionOfPseudoJet.hh"

using std::string;

struct Options {
  string name     = "job";
  string id       = "0";
  string input    = "";
  string out_dir  = "";
  string tow_list = "";
  string run_list = "";
  string triggers = "";
};

void ConvertTStarJetVector(TStarJetVectorContainer<TStarJetVector>* container, std::vector<fastjet::PseudoJet> & particles);

int main(int argc, char* argv[]) {
  
  // parse command line options
  Options opts;
  for (int i = 1; i < argc; ++i) {
    if(ParseStrFlag(string(argv[i]), "--out", &opts.out_dir) ||
       ParseStrFlag(string(argv[i]), "--input", &opts.input) ||
       ParseStrFlag(string(argv[i]), "--badRuns", &opts.run_list) ||
       ParseStrFlag(string(argv[i]), "--badTowers", &opts.tow_list) ||
       ParseStrFlag(string(argv[i]), "--triggers", &opts.triggers) ||
       ParseStrFlag(string(argv[i]), "--name", &opts.name) ||
       ParseStrFlag(string(argv[i]),   "--id", &opts.id))
      continue;
    std::cerr << "unknown argument: " << argv[i] << std::endl;
    return 1;
  }
  
  // build the input chain
  TChain* chain = new TChain("JetTree");
  if (HasEnding(opts.input, ".root"))
    chain->Add(opts.input.c_str());
  else if (HasEnding(opts.input, ".list"))
    chain = TStarJetPicoUtils::BuildChainFromFileList( opts.input.c_str() );
  else if (HasEnding(opts.input, ".txt"))
    chain = TStarJetPicoUtils::BuildChainFromFileList( opts.input.c_str() );
  else {
    std::cerr << "Unrecognized input file format: exiting." << std::endl;
    return 1;
  }
  
  // build output directory if it doesn't exist
  // using the boost::filesystem library
  boost::filesystem::path dir(opts.out_dir.c_str());
  boost::filesystem::create_directories(dir);
  string outfile_name = opts.out_dir + "/" + opts.name + opts.id + ".root";
  TFile out(outfile_name.c_str(), "RECREATE");
  
  // initialize the reader
  TStarJetPicoReader reader;
  reader.SetApplyFractionHadronicCorrection( true );
  reader.SetFractionHadronicCorrection( 0.999 );
  reader.SetRejectTowerElectrons( kFALSE );
  
  // set event cuts
  reader.SetInputChain(chain);
  reader.GetTowerCuts()->AddBadTowers(opts.tow_list.c_str());
  reader.GetTowerCuts()->SetMaxEtCut(1000);
  reader.GetTrackCuts()->SetDCACut(1.0);
  reader.GetTrackCuts()->SetMinNFitPointsCut(20);
  reader.GetTrackCuts()->SetFitOverMaxPointsCut(0.52);
  reader.GetTrackCuts()->SetMaxPtCut(1000);
  reader.GetEventCuts()->SetMaxEventPtCut(30);
  reader.GetEventCuts()->SetMaxEventEtCut(30);
  reader.GetEventCuts()->SetVertexZCut(30);
  reader.GetEventCuts()->SetTriggerSelection("HT");
  reader.GetEventCuts()->SetVertexZDiffCut(3);
  
  // get our bad run list
  if (opts.run_list != "") {
    std::ifstream bad_run_stream(opts.run_list);
    string value;
    while (bad_run_stream.good()) {
      getline ( bad_run_stream, value, ',' );
      string tmp = string( value, 0, value.length()-2 );
      tmp.erase(std::remove(tmp.begin(), tmp.end(), ' '), tmp.end());
      tmp.erase(std::remove(tmp.begin(), tmp.end(), ','), tmp.end());
      if (tmp.size() == 0) continue;
      reader.AddMaskedRun(stoi(tmp));
    }
  }
  
  reader.Init();
  
  // and the definition & selectors for fastjet
  const double resolution = 0.4;
  const double eta_acceptance_max   = 1.0;
  const double const_pt_min_high = 2.0;
  const double const_pt_min_low = 0.2;
  const double const_pt_max = 30.0;
  const double jet_eta_max = eta_acceptance_max - resolution;
  const double jet_pt_min = 10.0;
  const double jet_pt_max = 100.0;
  
  fastjet::JetDefinition jet_def = fastjet::JetDefinition(fastjet::antikt_algorithm, resolution);
  fastjet::Selector jet_selector    = fastjet::SelectorPtMin(jet_pt_min)   * fastjet::SelectorPtMax(jet_pt_max)   * fastjet::SelectorAbsRapMax(jet_eta_max);
  fastjet::Selector track_selector_high  = fastjet::SelectorPtMin(const_pt_min_high) * fastjet::SelectorPtMax(const_pt_max) * fastjet::SelectorAbsRapMax(eta_acceptance_max);
  fastjet::Selector track_selector_low  = fastjet::SelectorPtMin(const_pt_min_low) * fastjet::SelectorPtMax(const_pt_max) * fastjet::SelectorAbsRapMax(eta_acceptance_max);
  
  fastjet::JetDefinition bkg_jet_def = fastjet::JetDefinition(fastjet::kt_algorithm, resolution);
  fastjet::GhostedAreaSpec areaSpec( jet_eta_max + resolution, 1, 0.01);
  fastjet::AreaDefinition areaDef( fastjet::active_area_explicit_ghosts, areaSpec);
  
  fastjet::Selector bkg_selector = fastjet::SelectorAbsRapMax( eta_acceptance_max - resolution ) * (!fastjet::SelectorNHardest(2));
  
  // histograms
  TH2D* ajhigh = new TH2D("ajhigh", "ajhigh", 50, -0.6, 0.9, 800, -0.5, 799.5 );
  TH2D* ajlow = new TH2D("ajlow", "ajlow", 50, -0.6, 0.9, 800, -0.5, 799.5 );
  TH1D* ajhigh020 = new TH1D("ajhigh020", "ajhigh020", 50, -0.6, 0.9);
  TH1D* ajlow020 = new TH1D("ajlow020", "ajlow020", 50, -0.6, 0.9);
  TH1D* leadpt = new TH1D("leadpt", "leadpt", 100, 0, 100);
  TH1D* subpt = new TH1D("subpt", "subpt", 100, 0, 100);
  TH1D* leadptsoft = new TH1D("leadptsoft", "leadpt", 100, 0, 100);
  TH1D* subptsoft = new TH1D("subptsoft", "subpt", 100, 0, 100);
  TH1D* dR = new TH1D("dr", "dr", 100, 0, TMath::Pi()+0.2);
  
  while(reader.NextEvent()) {
    // Print out reader status every 10 seconds
    reader.PrintStatus(10);
    
    TStarJetPicoEvent* event = reader.GetEvent();
    TStarJetPicoEventHeader* header = event->GetHeader();
    
    // do jetfinding for hard core jets
    TStarJetVectorContainer<TStarJetVector>* container = reader.GetOutputContainer();
    std::vector<fastjet::PseudoJet> particles;
    ConvertTStarJetVector(container, particles);
    std::vector<fastjet::PseudoJet> high_const = track_selector_high(particles);
    std::vector<fastjet::PseudoJet> high_jets;
    std::vector<fastjet::PseudoJet> low_jets;
    if(high_const.size() != 0) {
      fastjet::ClusterSequence cluster(high_const, jet_def);
      high_jets = fastjet::sorted_by_pt(jet_selector(cluster.inclusive_jets()));
    }
    
    // if we have two that are matched back to back, pt >20, pt >10, match
    if (high_jets.size() < 2 ) continue;
    
    if (high_jets[0].pt() < 20 || fabs(fabs(high_jets[0].delta_phi_to(high_jets[1])) - TMath::Pi()) > 0.4) continue;
    
    // now recluster & do background subtraction
    std::vector<fastjet::PseudoJet> low_const = track_selector_low(particles);
    
    fastjet::ClusterSequenceArea cluster_low ( low_const, jet_def, areaDef );
    // Energy density estimate from median ( pt_i / area_i )
    fastjet::JetMedianBackgroundEstimator bkgdEstimator ( bkg_selector, bkg_jet_def, areaDef );
    bkgdEstimator.set_particles( low_const );
    // Subtract A*rho from the original pT
    fastjet::Subtractor bkgdSubtractor ( &bkgdEstimator );
    
    low_jets = fastjet::sorted_by_pt( bkgdSubtractor( cluster_low.inclusive_jets() ) );
    
    if (low_jets.size() < 2) continue;
    
    fastjet::Selector radial_selector = fastjet::SelectorCircle(resolution);
    radial_selector.set_reference(high_jets[0]);
    std::vector<fastjet::PseudoJet> matched_to_high = fastjet::sorted_by_pt(radial_selector(low_jets));
    
    radial_selector.set_reference(high_jets[1]);
    std::vector<fastjet::PseudoJet> matched_to_low = fastjet::sorted_by_pt(radial_selector(low_jets));
    
    if (matched_to_low.size() == 0 || matched_to_high.size() == 0) continue;
    
    // now we have our four jets: continue
    fastjet::PseudoJet lead_hard = high_jets[0];
    fastjet::PseudoJet sub_hard = high_jets[1];
    fastjet::PseudoJet lead_soft = matched_to_high[0];
    fastjet::PseudoJet sub_soft = matched_to_low[0];
    
    double aj_hard = (lead_hard.pt() - sub_hard.pt()) / (lead_hard.pt() + sub_hard.pt());
    double aj_soft = (lead_soft.pt() - sub_soft.pt()) / (lead_soft.pt() + sub_soft.pt());
    
    ajhigh->Fill(aj_hard, header->GetReferenceMultiplicity());
    ajlow->Fill(aj_soft, header->GetReferenceMultiplicity());
    dR->Fill(lead_hard.delta_phi_to(sub_hard));
    leadpt->Fill(lead_hard.pt());
    subpt->Fill(sub_hard.pt());
    leadptsoft->Fill(lead_soft.pt());
    subptsoft->Fill(sub_soft.pt());
    
    if (header->GetReferenceMultiplicity() > 269) {
      ajhigh020->Fill(aj_hard);
      ajlow020->Fill(aj_soft);
    }
    
  }
  
  
  out.Write();
  return 0;
}

void ConvertTStarJetVector(TStarJetVectorContainer<TStarJetVector>* container, std::vector<fastjet::PseudoJet> & particles) {
  
  // Transform TStarJetVectors into (FastJet) PseudoJets
  // ---------------------------------------------------
  TStarJetVector* sv;
  for(int i = 0; i < container->GetEntries() ; ++i) {
    sv = container->Get(i);
    
    fastjet::PseudoJet tmpPJ = fastjet::PseudoJet(*sv);
    tmpPJ.set_user_index(sv->GetCharge());
    particles.push_back(tmpPJ);
  }
}
