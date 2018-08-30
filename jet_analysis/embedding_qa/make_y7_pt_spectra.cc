#include "jet_analysis/util/arg_helper.hh"
#include "jet_analysis/util/trigger_lookup.hh"
#include "jet_analysis/util/reader_util.hh"
#include "jet_analysis/util/string_util.hh"
#include "jet_analysis/efficiency/run14_eff.hh"
#include "jet_analysis/efficiency/run7_eff.hh"
#include "jet_analysis/centrality/centrality_run14.hh"


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

// use glog and gflags
#include "gflags/gflags.h"
#include "glog/stl_logging.h"
#include "glog/logging.h"

DEFINE_string(name, "run7_pt", "output name");
DEFINE_string(input, "", "input data file/list");
DEFINE_string(outdir, "tmp", "output directory");
DEFINE_string(towList, "submit/y7_y6_bad_tower.txt", "bad tower list");
DEFINE_string(runList, "", "optional bad run list");
DEFINE_string(triggers, "y7mb", "trigger selection");

using std::string;

int main(int argc, char* argv[]) {
  
  string usage = "generate pT spectra for Run 7 data";
  
  gflags::SetUsageMessage(usage);
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);
 
  // check to make sure the input file exists
  if (!boost::filesystem::exists(FLAGS_input)) {
    std::cerr << "input file does not exist: " << FLAGS_input << std::endl;;
    return 1;
  }
  
  // build output directory if it doesn't exist, using boost::filesystem
  if (FLAGS_outdir.empty())
    FLAGS_outdir = "tmp";
  boost::filesystem::path dir(FLAGS_outdir.c_str());
  boost::filesystem::create_directories(dir);
  
  // build our input chain
  TChain* chain = NewChainFromInput(FLAGS_input);
  
  // initialize reader
  // create output file from the given directory, name & id
  string outfile_name = FLAGS_outdir + "/" + FLAGS_name + ".root";
  TFile out(outfile_name.c_str(), "RECREATE");
  
  // initialize the reader(s)
  TStarJetPicoReader* reader = new TStarJetPicoReader();
  InitReaderWithDefaults(reader, chain, FLAGS_towList, FLAGS_runList);
  
  // get the triggers IDs that will be used
  std::set<unsigned> triggers = GetTriggerIDs(FLAGS_triggers);
  
  Run7Eff* run7Eff = new Run7Eff();
  
  int cent_bins = 9;
  std::vector<std::pair<int, int>> CentBoundariesY7{{485, 1000}, {399, 484}, {269, 398}, {178, 268}, {114, 177},
                                                    {69, 113}, {39, 68}, {21, 38}, {10, 20}};
  
  // change to output file
  out.cd();
  
  // Histograms will calculate gaussian errors
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();
  
  // create pT histogram
  TH2D* nglobnprim = new TH2D("nglobnprim", "", 500, 0, 5000, 500, 0, 1500);
  TH2D* refmultnprim = new TH2D("refmultnprim", "", 800, 0, 800, 500, 0, 1500);
  TH2D* pt = new TH2D("pt", ";centrality;p_{T}", cent_bins, -0.5, cent_bins - 0.5, 50, 0, 5);
  TH2D* pt_corr = new TH2D("ptcorr", ";centrality;p_{T}",cent_bins, -0.5, cent_bins - 0.5, 50, 0, 5);
  TH2D* refmult = new TH2D("refmult", ";centrality;refmult", cent_bins, -0.5, cent_bins - 0.5, 800, 0, 800);
  TH2D* frac = new TH2D("discarded", "", cent_bins, -0.5, cent_bins - 0.5, 10, 0, 1.0);
  TH2D* nprim = new TH2D("nprim", "", cent_bins, -0.5, cent_bins - 0.5, 100, 0, 1200);
  TH2D* nsel = new TH2D("nsel", "", cent_bins, -0.5, cent_bins - 0.5, 100, 0, 12000);
  TH2D* nhitsfit = new TH2D("nhitsfit", "", cent_bins, -0.5, cent_bins - 0.5, 50, 0, 50);
  TH2D* nhitspos = new TH2D("nhitspos", "", cent_bins, -0.5, cent_bins - 0.5, 50, 0, 50);
  TH2D* nhitsfitfrac = new TH2D("nhitsfitfrac", "", cent_bins, -0.5, cent_bins - 0.5, 50, 0, 1.0);
  TH3D* dca = new TH3D("dcapt", "", cent_bins, -0.5, cent_bins - 0.5, 50, 0, 3, 50, 0, 5);
  
  std::vector<TProfile*> avg_eff(8);
  for (int i = 0; i < 8; ++i) {
    avg_eff[i] = new TProfile(MakeString("eff", i).c_str(), "", 100, 0, 5.0);
  }
  
  // start the event loop
  // --------------------
  while(reader->NextEvent()) {
    double counts = 0;
    double norm = 0;
    // Print out reader status every 10 seconds
    reader->PrintStatus(10);
    
    TStarJetPicoEvent* event = reader->GetEvent();
    TStarJetPicoEventHeader* header = event->GetHeader();
    
    nglobnprim->Fill(header->GetNGlobalTracks(), header->GetNOfPrimaryTracks());
    refmultnprim->Fill(header->GetReferenceMultiplicity(), header->GetNOfPrimaryTracks());
    
    int cent_bin = -1;
    // check if event fired a trigger we will use
    if (triggers.size() != 0) {
      bool use_event = false;
      for (auto trigger : triggers)
      if (header->HasTriggerId(trigger))
      use_event = true;
      if (!use_event) continue;
    }
    
    
    for (int i = 0; i < cent_bins; ++i) {
      if (header->GetGReferenceMultiplicity() >= CentBoundariesY7[i].first &&
          header->GetGReferenceMultiplicity() < CentBoundariesY7[i].second) {
        cent_bin = i;
        break;
      }
    }
    
    if (cent_bin < 0)
      continue;
    
    refmult->Fill(cent_bin, header->GetGReferenceMultiplicity());
    
    // get tracks & towers
    TList* tracks = reader->GetListOfSelectedTracks();
    int selected = 0;
    TIter nextTrack(tracks);
    while(TStarJetPicoPrimaryTrack* track = (TStarJetPicoPrimaryTrack*) nextTrack()) {
      
      if (fabs(track->GetEta()) > 1.0 || track->GetPt() < 0.2)
        continue;
      selected++;
      nhitsfit->Fill(cent_bin, track->GetNOfFittedHits());
      dca->Fill(cent_bin, track->GetDCA(), track->GetPt());
      nhitspos->Fill(cent_bin, track->GetNOfPossHits());
      nhitsfitfrac->Fill(cent_bin, (double)track->GetNOfFittedHits()/track->GetNOfPossHits());
      
      // do efficiency corrected pt spectrum
      double eff = -1;
      
      if (header->GetGReferenceMultiplicity() >= 485)
        eff = run7Eff->AuAuEff(track->GetPt(), track->GetEta(), 0);
      else if (header->GetGReferenceMultiplicity() >= 399)
        eff = run7Eff->AuAuEff(track->GetPt(), track->GetEta(), 1);
      else if (header->GetGReferenceMultiplicity() >= 269)
        eff = run7Eff->AuAuEff(track->GetPt(), track->GetEta(), 2);
      
      norm++;
      
      pt->Fill(cent_bin, track->GetPt());
      
      if (eff < 0 || eff > 1.0) {
        counts++;
        continue;
      }
      
      avg_eff[cent_bin]->Fill(track->GetPt(), eff);
      pt_corr->Fill(cent_bin, track->GetPt(), 1.0 / eff);
    }
    nsel->Fill(cent_bin, selected);
    nprim->Fill(cent_bin, selected);
    frac->Fill(cent_bin, counts/norm);
  }
  
  out.Write();
  out.Close();
  
  return 0;
}


