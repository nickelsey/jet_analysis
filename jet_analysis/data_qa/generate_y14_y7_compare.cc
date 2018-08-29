// extremely simple data comparison between y14 and y7
// to try and reduce any chances for errors......

#include "jet_analysis/util/arg_helper.hh"
#include "jet_analysis/util/trigger_lookup.hh"
#include "jet_analysis/util/reader_util.hh"
#include "jet_analysis/util/string_util.hh"
#include "jet_analysis/efficiency/run14_eff.hh"
#include "jet_analysis/efficiency/run7_eff.hh"
#include "jet_analysis/centrality/centrality_run14.hh"
#include "jet_analysis/util/root_print_routines.hh"

#include <string>
#include <iostream>

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TStyle.h"

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

using std::string;

DEFINE_string(y14Input, "/Users/nick/physics/data/y14/AuAu_200_MB_low_101_105_0.root", "y14 input file");
DEFINE_string(y7Input, "/Users/nick/physics/data/y7mb/newpicoDstcentralMB_8177020_DC4BA348C050D5562E7461357C4B341D_0.root", "y7 input file");
DEFINE_string(outdir, "y14_y7_pt_compare", "output directory");

DEFINE_string(y14TowList, "submit/y14_y6_bad_tower.txt", "bad tower list for y14");
DEFINE_string(y7TowList, "submit/y7_y6_bad_tower.txt", "bad tower list for y7");

DEFINE_string(y14RunList, "submit/y14_bad_run.txt", "bad run list for y14");
DEFINE_string(y7RunList, "", "bad run list for y7");

DEFINE_string(y14Triggers, "y14vpdmb30", "y14 trigger selection");
DEFINE_string(y7Triggers, "y7mb", "y7 trigger selection");

DEFINE_bool(central, true, "when true, only most central 0-20% events selected");

int main(int argc, char* argv[]) {
  
  string usage = "compares central pT spectra for y14 and y7";
  
  gflags::SetUsageMessage(usage);
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);
  
  // set drawing preferences for histograms and graphs
  gStyle->SetOptStat(false);
  gStyle->SetOptFit(false);
  gStyle->SetLegendBorderSize(0);
  
  // check to make sure the input file exists
  if (!boost::filesystem::exists(FLAGS_y14Input)) {
    std::cerr << "input file does not exist: " << FLAGS_y14Input << std::endl;;
    return 1;
  }
  if (!boost::filesystem::exists(FLAGS_y7Input)) {
    std::cerr << "input file does not exist: " << FLAGS_y7Input << std::endl;;
    return 1;
  }
  
  // build output directory if it doesn't exist, using boost::filesystem
  if (FLAGS_outdir.empty())
    FLAGS_outdir = "tmp";
  boost::filesystem::path dir(FLAGS_outdir.c_str());
  boost::filesystem::create_directories(dir);
  
  // make input chains
  TChain* y14_chain = NewChainFromInput(FLAGS_y14Input);
  TChain* y7_chain = NewChainFromInput(FLAGS_y7Input);
  
  // initialize the reader(s)
  TStarJetPicoReader* y14_reader = new TStarJetPicoReader();
  InitReaderWithDefaults(y14_reader, y14_chain, FLAGS_y14TowList, FLAGS_y14RunList);
  
  // tower cuts
  y14_reader->GetTowerCuts()->SetMaxEtCut(1000);            // essentially infinity - cut in eventcuts
  
  // track cuts
  y14_reader->GetTrackCuts()->SetDCACut(3.0);                // distance of closest approach to primary vtx
  y14_reader->GetTrackCuts()->SetMinNFitPointsCut(20);       // minimum fit points in track reco
  y14_reader->GetTrackCuts()->SetFitOverMaxPointsCut(0.0);  // minimum ratio of fit points used over possible
  y14_reader->GetTrackCuts()->SetMaxPtCut(1000);             // essentially infinity - cut in eventcuts
  
  // event cuts
  y14_reader->GetEventCuts()->SetMaxEventPtCut(30);          // Set Maximum track Pt
  y14_reader->GetEventCuts()->SetMaxEventEtCut(30);          // Set Maximum tower Et
  y14_reader->GetEventCuts()->SetVertexZCut(30);             // vertex z range (z = beam axis)
  y14_reader->GetEventCuts()->SetTriggerSelection("All");    // setting trigger selection - set to all, selected later
  y14_reader->GetEventCuts()->SetVertexZDiffCut(3);          // cut on Vz - VPD Vz
  
  
  TStarJetPicoReader* y7_reader = new TStarJetPicoReader();
  InitReaderWithDefaults(y7_reader, y7_chain, FLAGS_y7TowList, FLAGS_y7RunList);
  
  // tower cuts
  y7_reader->GetTowerCuts()->SetMaxEtCut(1000);            // essentially infinity - cut in eventcuts
  
  // track cuts
  y7_reader->GetTrackCuts()->SetDCACut(1.0);                // distance of closest approach to primary vtx
  y7_reader->GetTrackCuts()->SetMinNFitPointsCut(20);       // minimum fit points in track reco
  y7_reader->GetTrackCuts()->SetFitOverMaxPointsCut(0.52);  // minimum ratio of fit points used over possible
  y7_reader->GetTrackCuts()->SetMaxPtCut(1000);             // essentially infinity - cut in eventcuts
  
  // event cuts
  y7_reader->GetEventCuts()->SetMaxEventPtCut(30);          // Set Maximum track Pt
  y7_reader->GetEventCuts()->SetMaxEventEtCut(30);          // Set Maximum tower Et
  y7_reader->GetEventCuts()->SetVertexZCut(30);             // vertex z range (z = beam axis)
  y7_reader->GetEventCuts()->SetTriggerSelection("All");    // setting trigger selection - set to all, selected later
  y7_reader->GetEventCuts()->SetVertexZDiffCut(3);          // cut on Vz - VPD Vz
  
  // get triggers for y14 and y7
  std::set<unsigned> y14_triggers = GetTriggerIDs(FLAGS_y14Triggers);
  std::set<unsigned> y7_triggers = GetTriggerIDs(FLAGS_y7Triggers);
  
  // initialize centrality definition for run 14
  // and set a lower limit on 0-20% for year 7
  CentralityRun14 centrality;
  
  unsigned y7_grefmult_lower_bound = 269;
  
  // Histograms will calculate gaussian errors
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();
  
  // define histograms
  TH1D* y7_pt = new TH1D("y7pt", ";p_{T}", 20, 0, 10);
  TH1D* y14_pt = new TH1D("y14pt", ";p_{T}", 20, 0, 10);
  
  TH1D* y7_refmult = new TH1D("y7refmult", ";refmult", 100, 0, 800);
  TH1D* y14_refmult = new TH1D("y14refmult", ";refmult", 100, 0, 800);
  
  TH1D* y7_nprim = new TH1D("y7nprim", ";N_{primary}", 100, 0, 1000);
  TH1D* y14_nprim = new TH1D("y14nprim", ";N_{primary}", 100, 0, 1000);
  
  TH2D* y14_rank_dca = new TH2D("y14rankdca", ";rank;dca", 13, -10, 3, 50, 0, 3);
  TH2D* y14_rank_nhit = new TH2D("y14ranknhit", ";rank;nhit", 13, -10, 3, 50, 0, 50);
  TH2D* y14_rank_nhitpos = new TH2D("y14ranknhitpos", ";rank;nhitpos", 13, -10, 3, 50, 0, 50);
  TH2D* y14_rank_nprim = new TH2D("y14ranknprim", ";rank;nprim", 13, -10, 3, 100, 0, 1000);
  
  TH2D* y14_nglob_dca = new TH2D("y14nglobdca", ";N_{global};dca", 40, 0, 4000, 50, 0, 3);
  TH2D* y14_nglob_nhit = new TH2D("y14nglobnhit", ";N_{global};nhit", 40, 0, 4000, 50, 0, 50);
  TH2D* y14_nglob_nhitpos = new TH2D("y14nglobnhitpos", ";N_{global};nhitpos", 40, 0, 4000, 50, 0, 50);
  TH2D* y14_nglob_nprim = new TH2D("y14nglobnprim", ";N_{global};nprim", 40, 0, 4000, 100, 0, 1000);
  
  while(y14_reader->NextEvent()) {
    // Print out reader status every 10 seconds
    y14_reader->PrintStatus(10);
    
    TStarJetPicoEvent* event = y14_reader->GetEvent();
    TStarJetPicoEventHeader* header = event->GetHeader();
    
    if (y14_triggers.size() != 0) {
      bool use_event = false;
      for (auto trigger : y14_triggers)
        if (header->HasTriggerId(trigger))
          use_event = true;
      if (!use_event) continue;
    }
    
    centrality.setEvent(header->GetRunId(), header->GetReferenceMultiplicity(),
                        header->GetZdcCoincidenceRate(), header->GetPrimaryVertexZ());
    if (centrality.centrality16() < 0 || centrality.centrality16() > 15)
      continue;
    
    if (FLAGS_central && centrality.centrality16() > 3)
      continue;
    
    if (FLAGS_central && header->GetGReferenceMultiplicity() < y7_grefmult_lower_bound)
      continue;
    
    y14_refmult->Fill(header->GetReferenceMultiplicity());
    
    // get tracks & towers
    TList* tracks = y14_reader->GetListOfSelectedTracks();
    TIter nextTrack(tracks);
    int nprim = 0;
    while(TStarJetPicoPrimaryTrack* track = (TStarJetPicoPrimaryTrack*) nextTrack()) {
      if (fabs(track->GetEta()) > 1.0 || track->GetPt() < 0.2)
        continue;
      
      nprim++;
      
      y14_pt->Fill(track->GetPt());
      
      y14_rank_dca->Fill(header->GetPrimaryVertexRanking(), track->GetDCA());
      y14_rank_nhit->Fill(header->GetPrimaryVertexRanking(), track->GetNOfFittedHits());
      y14_rank_nhitpos->Fill(header->GetPrimaryVertexRanking(), track->GetNOfPossHits());
      
      y14_nglob_dca->Fill(header->GetNGlobalTracks(), track->GetDCA());
      y14_nglob_nhit->Fill(header->GetNGlobalTracks(), track->GetNOfFittedHits());
      y14_nglob_nhitpos->Fill(header->GetNGlobalTracks(), track->GetNOfPossHits());
      
    }
    y14_nprim->Fill(nprim);
    y14_rank_nprim->Fill(header->GetPrimaryVertexRanking(), nprim);
    y14_nglob_nprim->Fill(header->GetNGlobalTracks(), nprim);
  }
  
  while (y7_reader->NextEvent()) {
    // Print out reader status every 10 seconds
    y7_reader->PrintStatus(10);
    
    TStarJetPicoEvent* event = y7_reader->GetEvent();
    TStarJetPicoEventHeader* header = event->GetHeader();
    
    if (y7_triggers.size() != 0) {
      bool use_event = false;
      for (auto trigger : y7_triggers)
        if (header->HasTriggerId(trigger))
          use_event = true;
      if (!use_event) continue;
    }
    
    if (FLAGS_central && header->GetGReferenceMultiplicity() < y7_grefmult_lower_bound)
      continue;
    
    y7_refmult->Fill(header->GetGReferenceMultiplicity());
    
    TList* tracks = y7_reader->GetListOfSelectedTracks();
    TIter nextTrack(tracks);
    int nprim = 0;
    while(TStarJetPicoPrimaryTrack* track = (TStarJetPicoPrimaryTrack*) nextTrack()) {
      if (fabs(track->GetEta()) > 1.0 || track->GetPt() < 0.2)
        continue;
      
      nprim++;
      
      y7_pt->Fill(track->GetPt());
    }
    y7_nprim->Fill(nprim);
  }
  
  // now normalize
  y7_pt->Scale(1.0 / y7_nprim->GetEntries());
  y14_pt->Scale(1.0 / y14_nprim->GetEntries());
  
  y7_refmult->Scale(1.0 / y7_refmult->Integral());
  y14_refmult->Scale(1.0 / y14_refmult->Integral());
  
  y7_nprim->Scale(1.0 / y7_nprim->Integral());
  y14_nprim->Scale(1.0 / y14_nprim->Integral());
  
  histogramOpts hopts;
  canvasOpts copts;
  copts.do_legend = false;
  canvasOpts coptslogy;
  coptslogy.do_legend = false;
  coptslogy.log_y = true;
  
  Overlay1D(y7_refmult, y14_refmult, "run 7", "run 14", hopts, coptslogy, FLAGS_outdir, "refmult",
            "", "refmult", "fraction", "");
  Overlay1D(y7_nprim, y14_nprim, "run 7", "run 14", hopts, coptslogy, FLAGS_outdir, "nprim",
            "", "N_{primary}", "fraction", "");
  PrintWithRatio(y7_pt, y14_pt, "run 7", "run 14", hopts, coptslogy, FLAGS_outdir, "pt", "", "p_{T}",
                 "1/N_{events}dN/dp_{T}", "");
  
  TProfile* prof_dca = y14_rank_dca->ProfileX();
  TH1D* proj_dca = y14_rank_dca->ProjectionX();
  TProfile* prof_nhit = y14_rank_nhit->ProfileX();
  TH1D* proj_nhit = y14_rank_nhit->ProjectionX();
  TProfile* prof_nhitposs = y14_rank_nhitpos->ProfileX();
  TH1D* proj_nhitposs = y14_rank_nhitpos->ProjectionX();
  TProfile* prof_nprim = y14_rank_nprim->ProfileX();
  TH1D* proj_nprim = y14_rank_nprim->ProjectionX();
  
  TProfile* nglob_prof_dca = y14_nglob_dca->ProfileX();
  TH1D* nglob_proj_dca = y14_nglob_dca->ProjectionX();
  TProfile* nglob_prof_nhit = y14_nglob_nhit->ProfileX();
  TH1D* nglob_proj_nhit = y14_nglob_nhit->ProjectionX();
  TProfile* nglob_prof_nhitposs = y14_nglob_nhitpos->ProfileX();
  TH1D* nglob_proj_nhitposs = y14_nglob_nhitpos->ProjectionX();
  TProfile* nglob_prof_nprim = y14_nglob_nprim->ProfileX();
  TH1D* nglob_proj_nprim = y14_nglob_nprim->ProjectionX();
  
  PrettyPrint1D(prof_dca, hopts, copts, "", FLAGS_outdir, "dca", "", "rank", "<DCA>");
  PrettyPrint1D(prof_nhit, hopts, copts, "", FLAGS_outdir, "nhit", "", "rank", "<N_{hit}>");
  PrettyPrint1D(prof_nhitposs, hopts, copts, "", FLAGS_outdir, "nhitposs", "", "rank", "<N_{hitposs}>");
  PrettyPrint1D(prof_nprim, hopts, copts, "", FLAGS_outdir, "nprimprof", "", "rank", "<N_{primary}");
  
  PrettyPrint1D(proj_dca, hopts, copts, "", FLAGS_outdir, "dcaproj", "", "rank", "count");
  PrettyPrint1D(proj_nhit, hopts, copts, "", FLAGS_outdir, "nhitproj", "", "rank", "count");
  PrettyPrint1D(proj_nhitposs, hopts, copts, "", FLAGS_outdir, "nhitpossproj", "", "rank", "count");
  PrettyPrint1D(proj_nprim, hopts, copts, "", FLAGS_outdir, "nprimproj", "", "rank", "count");
  
  PrettyPrint1D(nglob_prof_dca, hopts, copts, "", FLAGS_outdir, "nglobdca", "", "N_{global}", "<DCA>");
  PrettyPrint1D(nglob_prof_nhit, hopts, copts, "", FLAGS_outdir, "nglobnhit", "", "N_{global}", "<N_{hit}>");
  PrettyPrint1D(nglob_prof_nhitposs, hopts, copts, "", FLAGS_outdir, "nglobnhitposs", "", "N_{global}", "<N_{hitposs}>");
  PrettyPrint1D(nglob_prof_nprim, hopts, copts, "", FLAGS_outdir, "nglobnprimprof", "", "N_{global}", "<N_{primary}");
  
  PrettyPrint1D(nglob_proj_dca, hopts, copts, "", FLAGS_outdir, "nglobdcaproj", "", "N_{global}", "count");
  PrettyPrint1D(nglob_proj_nhit, hopts, copts, "", FLAGS_outdir, "nglobnhitproj", "", "N_{global}", "count");
  PrettyPrint1D(nglob_proj_nhitposs, hopts, copts, "", FLAGS_outdir, "nglobnhitpossproj", "", "N_{global}", "count");
  PrettyPrint1D(nglob_proj_nprim, hopts, copts, "", FLAGS_outdir, "nglobnprimproj", "", "N_{global}", "count");
  
 
  return 0;
}
