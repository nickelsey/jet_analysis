#include "jet_analysis/util/arg_helper.hh"
#include "jet_analysis/util/trigger_lookup.hh"
#include "jet_analysis/util/reader_util.hh"
#include "jet_analysis/util/string_util.hh"
#include "jet_analysis/efficiency/run14_eff.hh"
#include "jet_analysis/efficiency/run7_eff.hh"
#include "jet_analysis/centrality/centrality_run14.hh"


#include <string>

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

using std::string;

struct Options {
  string name        = "job";    /* output file name */
  string id          = "0";      /* job id */
  string input       = "";       /* root file/root file list*/
  string out_dir     = "tmp";    /* directory to save output in */
  string tow_list    = "";       /* list of hot towers to remove */
  string run_list    = "";       /* list of runs to remove */
  string triggers    = "";       /* triggers to analyze (see trigger_lookup.hh) */
  string effcurves   = "";       /* file holding y14 efficiency curves */
  bool useY14Eff     = false;    /* use y14 efficiency curves when doing efficiency corrections */
  bool useY7Eff      = false;    /* use y7 efficiency curves when doing efficiency corrections */
  int nhitsfit       = 10;       /* nhits fit cut */
  double eta         = 1.0;      /* eta cut for tracks */
  double dca         = 3.0;      /* track dca cut */
  double fitfrac     = 0.0;      /* nhits fit/ nhits prossible */
};

int main(int argc, char* argv[]) {
  
  // parse command line options
  // --------------------------
  Options opts;
  for (int i = 1; i < argc; ++i) {
    if (ParseStrFlag(string(argv[i]), "--name", &opts.name) ||
        ParseStrFlag(string(argv[i]), "--id", &opts.id) ||
        ParseStrFlag(string(argv[i]), "--input", &opts.input) ||
        ParseStrFlag(string(argv[i]), "--outDir", &opts.out_dir) ||
        ParseStrFlag(string(argv[i]), "--towList", &opts.tow_list) ||
        ParseStrFlag(string(argv[i]), "--runList", &opts.run_list) ||
        ParseStrFlag(string(argv[i]), "--triggers", &opts.triggers) ||
        ParseStrFlag(string(argv[i]), "--efficiencyCurves", &opts.effcurves) ||
        ParseBoolFlag(string(argv[i]), "--year7", &opts.useY7Eff) ||
        ParseBoolFlag(string(argv[i]), "--year14", &opts.useY14Eff) ||
        ParseIntFlag(string(argv[i]), "--nhitsfit", &opts.nhitsfit) ||
        ParseFloatFlag(string(argv[i]), "--eta", &opts.eta) ||
        ParseFloatFlag(string(argv[i]), "--dca", &opts.dca) ||
        ParseFloatFlag(string(argv[i]), "--fitfrac", &opts.fitfrac)) continue;
    std::cerr << "Unknown command line option: " << argv[i] << std::endl;
    return 1;
  }
 
  if (opts.useY7Eff == false && opts.useY14Eff == false) {
    std::cerr << "need to select year 14 or year 7" << std::endl;
    return 1;
  }
  
  // check to make sure the input file exists
  if (!boost::filesystem::exists(opts.input)) {
    std::cerr << "input file does not exist: " << opts.input << std::endl;;
    return 1;
  }
  
  // build output directory if it doesn't exist, using boost::filesystem
  if (opts.out_dir.empty())
    opts.out_dir = "tmp";
  boost::filesystem::path dir(opts.out_dir.c_str());
  boost::filesystem::create_directories(dir);
  
  // build our input chain
  TChain* chain = NewChainFromInput(opts.input);
  
  // create output file from the given directory, name & id
  string outfile_name = opts.out_dir + "/" + opts.name + opts.id + ".root";
  TFile out(outfile_name.c_str(), "RECREATE");
  
  // initialize the reader(s)
  TStarJetPicoReader* reader = new TStarJetPicoReader();
  InitReaderWithDefaults(reader, chain, opts.tow_list, opts.run_list);
  reader->GetTrackCuts()->SetDCACut(opts.dca);                // distance of closest approach to primary vtx
  reader->GetTrackCuts()->SetMinNFitPointsCut(opts.nhitsfit);       // minimum fit points in track reco
  reader->GetTrackCuts()->SetFitOverMaxPointsCut(opts.fitfrac);  // minimum ratio of fit points used over possible
  reader->GetTrackCuts()->SetMaxPtCut(1000);             // essentially infinity - cut in eventcuts
  
  
  // get the triggers IDs that will be used
  std::set<unsigned> triggers = GetTriggerIDs(opts.triggers);
  
  Run14Eff* run14Eff;
  if (opts.useY14Eff) {
    run14Eff = new Run14Eff();
    run14Eff->loadFile(opts.effcurves);
  }
  
  Run7Eff* run7Eff;
  if (opts.useY7Eff) {
    run7Eff = new Run7Eff();
    run7Eff->loadFile(opts.effcurves);
  }
  
  // initialize centrality definition for run 14
  // and set a lower limit on 0-20% for year 7
  CentralityRun14 centrality;
  //int year7_centrality_cut = 269;
  int year7_centrality_cut = 485;
  // define upper and lower cut for 40-50%
  //int year7_centrality_cut_upper = 114;
  //int year7_centrality_cut_lower = 69;
  
  // change to output file
  out.cd();
  
  // Histograms will calculate gaussian errors
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();
  
  // create pT histogram
  TH1D* pt = new TH1D("pt", ";p_{T}", 100, 0, 5);
  TH1D* pt_corr = new TH1D("ptcorr", ";p_{T}", 100, 0, 5);
  TH1D* refmult = new TH1D("refmult", ";refmult", 800, 0, 800);
  TH1D* frac = new TH1D("discarded", "", 10, 0, 1.0);
  TH1D* nprim = new TH1D("nprim", "", 100, 0, 2000);
  TH1D* nsel = new TH1D("nsel", "", 100, 0, 2000);
  TProfile* avg_eff = new TProfile("eff", "", 100, 0, 5.0);
  TH1D* nhitsfit = new TH1D("nhitsfit", "", 50, 0, 50);
  TH1D* nhitspos = new TH1D("nhitspos", "", 50, 0, 50);
  TH1D* nhitsfitfrac = new TH1D("nhitsfitfrac", "", 50, 0, 1.0);
  TH2D* dca = new TH2D("dcapt", "", 50, 0, 3, 50, 0, 5);
  
  TH1D* nhitpos_outer_vz = new TH1D("outervz", "", 60, 0, 60);
  TH1D* nhitpos_inner_vz = new TH1D("innervz", "", 60, 0, 60);
  
  // 2D eta/phi histograms
  TH2D* etaphi_pt0 = new TH2D("etaphi0", ";#eta;#phi", 40, -1, 1, 40, -TMath::Pi(), TMath::Pi());
  TH2D* etaphi_pt1 = new TH2D("etaphi1", ";#eta;#phi", 40, -1, 1, 40, -TMath::Pi(), TMath::Pi());
  TH2D* etaphi_pt2 = new TH2D("etaphi2", ";#eta;#phi", 40, -1, 1, 40, -TMath::Pi(), TMath::Pi());
  TH2D* etaphi_pt3 = new TH2D("etaphi3", ";#eta;#phi", 40, -1, 1, 40, -TMath::Pi(), TMath::Pi());
  
  // start the event loop
  // --------------------
  while(reader->NextEvent()) {
    double counts = 0;
    double norm = 0;
    // Print out reader status every 10 seconds
    reader->PrintStatus(10);
    
    TStarJetPicoEvent* event = reader->GetEvent();
    TStarJetPicoEventHeader* header = event->GetHeader();
    
    
    int cent_bin = 0;
    // check if event fired a trigger we will use
    if (triggers.size() != 0) {
      bool use_event = false;
      for (auto trigger : triggers)
        if (header->HasTriggerId(trigger))
          use_event = true;
      if (!use_event) continue;
    }
    
    if (opts.useY7Eff) {
      if (header->GetGReferenceMultiplicity() < year7_centrality_cut)
        continue;
      refmult->Fill(header->GetGReferenceMultiplicity());
    }
    
    if (opts.useY14Eff) {
      centrality.setEvent(header->GetRunId(), header->GetReferenceMultiplicity(),
                          header->GetZdcCoincidenceRate(), header->GetPrimaryVertexZ());
      if (centrality.centrality16() > 0 || centrality.centrality16() < 0)
        continue;
      refmult->Fill(header->GetReferenceMultiplicity());
      cent_bin = centrality.centrality16();
    }
    
    TStarJetVectorContainer<TStarJetVector>* container = reader->GetOutputContainer();
    
    // get tracks & towers
    TList* tracks = reader->GetListOfSelectedTracks();
    int selected = 0;
    TIter nextTrack(tracks);
    while(TStarJetPicoPrimaryTrack* track = (TStarJetPicoPrimaryTrack*) nextTrack()) {
  
      if (fabs(track->GetEta()) > opts.eta || track->GetPt() < 0.2)
        continue;
      selected++;
      nhitsfit->Fill(track->GetNOfFittedHits());
      dca->Fill(track->GetDCA(), track->GetPt());
      nhitspos->Fill(track->GetNOfPossHits());
      nhitsfitfrac->Fill((double)track->GetNOfFittedHits()/track->GetNOfPossHits());
      
      if (fabs(header->GetPrimaryVertexZ()) < 10)
        nhitpos_inner_vz->Fill(track->GetNOfPossHits());
      if (fabs(header->GetPrimaryVertexZ()) > 25)
        nhitpos_outer_vz->Fill(track->GetNOfPossHits());
      
      if (track->GetPt() > 0.2 && track->GetPt() < 0.5)
        etaphi_pt0->Fill(track->GetEta(), track->GetPhi());
      else if (track->GetPt() < 1.0)
        etaphi_pt1->Fill(track->GetEta(), track->GetPhi());
      else if (track->GetPt() < 3.0)
        etaphi_pt2->Fill(track->GetEta(), track->GetPhi());
      else {
        etaphi_pt3->Fill(track->GetEta(), track->GetPhi());
      }
    }
    nsel->Fill(selected);
    
    TStarJetVector* sv;
    if (opts.useY7Eff) {
      for (int i = 0; i < container->GetEntries(); ++i) {
        sv = container->Get(i);
        if (sv->GetCharge() && fabs(sv->Eta()) < opts.eta && sv->Pt() > 0.2) {
          double eff = run7Eff->AuAuEff(sv->Pt(), sv->Eta(), 0);
          norm++;
          if (eff <= 0.0 || eff > 1.0) {
            counts++;
            continue;
          }
          pt->Fill(sv->Pt());
          pt_corr->Fill(sv->Pt(), 1.0 / eff);
          avg_eff->Fill(sv->Pt(), eff);
          
        }
      }
      nprim->Fill(norm);
      frac->Fill(counts/norm);
      
    }
    else if (opts.useY14Eff) {
      for (int i = 0; i < container->GetEntries(); ++i) {
        sv = container->Get(i);
        if (sv->GetCharge() && fabs(sv->Eta()) < opts.eta && sv->Pt() > 0.2) {
          double eff = run14Eff->AuAuEff(sv->Pt(), sv->Eta(), cent_bin, header->GetZdcCoincidenceRate());
          norm++;
          if (eff <= 0.0 || eff > 1.0) {
            counts++;
            continue;
          }
          pt->Fill(sv->Pt());
          pt_corr->Fill(sv->Pt(), 1.0 / eff);
          avg_eff->Fill(sv->Pt(), eff);
        }
      }
      nprim->Fill(norm);
      frac->Fill(counts/norm);
    }
  }
  
  out.Write();
  out.Close();
  return 0;
}
