#include "jet_analysis/util/trigger_lookup.hh"
#include "jet_analysis/util/reader_util.hh"
#include "jet_analysis/util/string_util.hh"
#include "jet_analysis/centrality/centrality_run14.hh"
#include "jet_analysis/util/common.hh"
#include "jet_analysis/util/vector_conversion.hh"

#include <string>
#include <iostream>

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"

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

DEFINE_string(input, "/Users/nick/physics/data/y14/AuAu_200_MB_low_101_105_0.root", "y14 input file");
DEFINE_string(outName, "bkg_rho", "Name for output root file");
DEFINE_string(outdir, "min_bias_rho", "output directory");
DEFINE_string(towlist, "submit/y14_y6_bad_tower.txt", "bad tower list");

int main(int argc, char* argv[]) {
  
  string usage = "Calculates background rho for 0-10%, for Dan";
  
  gflags::SetUsageMessage(usage);
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);

   // build output directory if it doesn't exist, using boost::filesystem
  if (FLAGS_outdir.empty()) FLAGS_outdir = "tmp";
  boost::filesystem::path dir(FLAGS_outdir.c_str());
  boost::filesystem::create_directories(dir);

  if (!boost::filesystem::exists(FLAGS_input)) {
    std::cerr << "input file does not exist: " << FLAGS_input << std::endl;;
    return 1;
  }

  // check to see which dataset is being used
  bool is_y14 = false;
  bool is_y11 = false;
  bool is_y7  = false;
  if (FLAGS_input.find("y11") != std::string::npos)
    is_y11 = true;
  if (FLAGS_input.find("y14") != std::string::npos)
    is_y14 = true;
  if (FLAGS_input.find("y7") != std::string::npos)
    is_y7 = true;

  if (!is_y14 && !is_y7 && !is_y11) {
    LOG(ERROR) << "Can't identify input dataset";
    return 1;
  }

  if ((is_y14 && is_y7) || (is_y14 && is_y11) || (is_y11 && is_y7)) {
    LOG(ERROR) << "input is ambiguous, can't distinguish between y7, y11, y14";
    return 1;
  }

  // build our input chain
  TChain* chain = NewChainFromInput(FLAGS_input);

  // initialize the reader(s)
  TStarJetPicoReader* reader = new TStarJetPicoReader();
  InitReaderWithDefaults(reader, chain, FLAGS_towlist, "");
  reader->GetTrackCuts()->SetDCACut(3.0);                // distance of closest approach to primary vtx
  reader->GetTrackCuts()->SetMinNFitPointsCut(20);       // minimum fit points in track reco
  reader->GetTrackCuts()->SetFitOverMaxPointsCut(0.52);  // minimum ratio of fit points used over possible
  reader->GetTrackCuts()->SetMaxPtCut(1000);             // essentially infinity - cut in eventcuts
  
  
	// create output file from the given directory, name & id
  string outfile_name = FLAGS_outdir + "/" + FLAGS_outName + ".root";
  TFile out(outfile_name.c_str(), "RECREATE");


  TH1D* rho = new TH1D("rho", "", 100, 0, 100);
  TH1D* refmult = new TH1D("ref", "", 100, 0, 800);

  // get the triggers IDs that will be used
  std::set<unsigned> triggers;
  if (is_y14)
    triggers = GetTriggerIDs("y14vpdmb30");
  else if (is_y7)
    triggers = GetTriggerIDs("y7mb");
  else if (is_y11)
    triggers = GetTriggerIDs("y11mb");

  // create centrality definition
  CentralityRun14 centrality;

  // fastjet stuff
  fastjet::Selector track_pt_min_selector = fastjet::SelectorAbsEtaMax(1.0) * fastjet::SelectorPtMin(0.2);
  fastjet::Selector jet_selector = fastjet::SelectorAbsEtaMax(0.6);
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 0.4);

  fastjet::JetDefinition bkg_jet_def(fastjet::kt_algorithm, 0.4);
  fastjet::Selector bkg_selector = fastjet::SelectorAbsRapMax(0.6);  // * (!fastjet::SelectorNHardest(2));
  fastjet::GhostedAreaSpec areaSpec( 1.4, 1, 0.01);
  fastjet::AreaDefinition areaDef( fastjet::active_area_explicit_ghosts, areaSpec);
  // start the event loop
  // --------------------
  while(reader->NextEvent()) {
    // Print out reader status every 10 seconds
    reader->PrintStatus(10);
    
    TStarJetPicoEvent* event = reader->GetEvent();
    TStarJetPicoEventHeader* header = event->GetHeader();
    
    // check if event fired a trigger we will use
    if (triggers.size() != 0) {
      bool use_event = false;
      for (auto trigger : triggers)
        if (header->HasTriggerId(trigger))
          use_event = true;
      if (!use_event) continue;
    }
    
    int cent = -1;

    if (is_y14) {
      centrality.setEvent(header->GetRunId(), header->GetReferenceMultiplicity(),
                        header->GetZdcCoincidenceRate(), header->GetPrimaryVertexZ());
      cent = centrality.centrality9();
    }
    else if (is_y11) {
      cent = abs(8 - header->GetReferenceCentrality());
    }
    else if (is_y7) {
      if (header->GetProperReferenceMultiplicity() >= 399)
        cent = 0;
    }
    
    if (cent < 0 || cent >= 2)
      continue;
    
      // get the vector container
      TStarJetVectorContainer<TStarJetVector>* container = reader->GetOutputContainer();
      std::vector<fastjet::PseudoJet> primary_particles;
      ConvertTStarJetVector(container, primary_particles);
      
      // select tracks above the minimum pt threshold
      primary_particles = track_pt_min_selector(primary_particles);
    
    // Energy density estimate from median ( pt_i / area_i )
    fastjet::JetMedianBackgroundEstimator bkgdEstimator ( bkg_selector, bkg_jet_def, areaDef );
    bkgdEstimator.set_particles(primary_particles);

    rho->Fill(bkgdEstimator.rho());
    refmult->Fill(header->GetProperReferenceMultiplicity());

  }

  rho->Write();
  refmult->Write();
  out.Close();

  return 0;
}