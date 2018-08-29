#include "jet_analysis/util/arg_helper.hh"
#include "jet_analysis/util/trigger_lookup.hh"
#include "jet_analysis/util/reader_util.hh"
#include "jet_analysis/util/string_util.hh"
#include "jet_analysis/util/root_print_routines.hh"
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


DEFINE_string(name, "track_quality", "output name");
DEFINE_string(data, "/Users/nick/physics/data/p18if/AuAu_200_MB_low_110_0.root", "input file");
DEFINE_string(outdir, "trackquality", "output directory");
DEFINE_string(TowList, "submit/y14_y6_bad_tower.txt", "bad tower list");
DEFINE_string(RunList, "submit/y14_bad_run.txt", "bad run list");
DEFINE_string(triggers, "y14vpdmb30", "trigger selection");

DEFINE_double(ptmin, 0.2, "min pt");

int main(int argc, char* argv[]) {
  string usage = "compares productions of run 14 on an event-by-event basis";
  
  gflags::SetUsageMessage(usage);
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);
  
  // set drawing preferences for histograms and graphs
  gStyle->SetOptStat(false);
  gStyle->SetOptFit(false);
  gStyle->SetLegendBorderSize(0);
  
  if (!boost::filesystem::exists(FLAGS_data)) {
    std::cerr << "input file does not exist: " << FLAGS_data << std::endl;;
    return 1;
  }
  
  // build output directory if it doesn't exist, using boost::filesystem
  if (FLAGS_outdir.empty())
    FLAGS_outdir = "tmp";
  boost::filesystem::path dir(FLAGS_outdir.c_str());
  boost::filesystem::create_directories(dir);
  
  // build our input chain
  TChain* chain = NewChainFromInput(FLAGS_data);
  
  // create output file from the given directory, name & id
  string outfile_name = FLAGS_outdir + "/" + FLAGS_name + ".root";
  TFile out(outfile_name.c_str(), "RECREATE");
  
  // initialize the reader(s)
  TStarJetPicoReader* reader = new TStarJetPicoReader();
  InitReaderWithDefaults(reader, chain, FLAGS_TowList, FLAGS_RunList);
  
  // tower cuts
  reader->GetTowerCuts()->SetMaxEtCut(1000);            // essentially infinity - cut in eventcuts
  
  // track cuts
  reader->GetTrackCuts()->SetDCACut(1.0);                // distance of closest approach to primary vtx
  reader->GetTrackCuts()->SetMinNFitPointsCut(15);       // minimum fit points in track reco
  reader->GetTrackCuts()->SetFitOverMaxPointsCut(0.0);  // minimum ratio of fit points used over possible
  reader->GetTrackCuts()->SetMaxPtCut(1000);             // essentially infinity - cut in eventcuts
  
  // event cuts
  reader->GetEventCuts()->SetMaxEventPtCut(30);          // Set Maximum track Pt
  reader->GetEventCuts()->SetMaxEventEtCut(30);          // Set Maximum tower Et
  reader->GetEventCuts()->SetVertexZCut(30);             // vertex z range (z = beam axis)
  reader->GetEventCuts()->SetTriggerSelection("All");    // setting trigger selection - set to all, selected later
  reader->GetEventCuts()->SetVertexZDiffCut(3);          // cut on Vz - VPD Vz
  
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();
  
  // histograms
  std::vector<TH1D*> pt_pos_near;
  std::vector<TH1D*> nhit_pos_near;
  std::vector<TH1D*> dca_pos_near;
  std::vector<TH1D*> pt_neg_near;
  std::vector<TH1D*> nhit_neg_near;
  std::vector<TH1D*> dca_neg_near;
  
  std::vector<TH1D*> pt_pos_far;
  std::vector<TH1D*> nhit_pos_far;
  std::vector<TH1D*> dca_pos_far;
  std::vector<TH1D*> pt_neg_far;
  std::vector<TH1D*> nhit_neg_far;
  std::vector<TH1D*> dca_neg_far;
  
  for (int i = 0; i < 9; ++i) {
    pt_pos_near.push_back(new TH1D(MakeString("pt", i, "posnear").c_str(), ";p_{T}", 15, 0, 5));
    nhit_pos_near.push_back(new TH1D(MakeString("nhit", i, "posnear").c_str(), ";p_{T}", 50, 0, 50));
    dca_pos_near.push_back(new TH1D(MakeString("dca", i, "posnear").c_str(), ";p_{T}", 20, 0, 2));
    pt_neg_near.push_back(new TH1D(MakeString("pt", i, "negnear").c_str(), ";p_{T}", 15, 0, 5));
    nhit_neg_near.push_back(new TH1D(MakeString("nhit", i, "negnear").c_str(), ";p_{T}", 50, 0, 50));
    dca_neg_near.push_back(new TH1D(MakeString("dca", i, "negnear").c_str(), ";p_{T}", 20, 0, 2));
    
    pt_pos_far.push_back(new TH1D(MakeString("pt", i, "posfar").c_str(), ";p_{T}", 15, 0, 5));
    nhit_pos_far.push_back(new TH1D(MakeString("nhit", i, "posfar").c_str(), ";p_{T}", 50, 0, 50));
    dca_pos_far.push_back(new TH1D(MakeString("dca", i, "posfar").c_str(), ";p_{T}", 20, 0, 2));
    pt_neg_far.push_back(new TH1D(MakeString("pt", i, "negfar").c_str(), ";p_{T}", 15, 0, 5));
    nhit_neg_far.push_back(new TH1D(MakeString("nhit", i, "negfar").c_str(), ";p_{T}", 50, 0, 50));
    dca_neg_far.push_back(new TH1D(MakeString("dca", i, "negfar").c_str(), ";p_{T}", 20, 0, 2));
    
  }
  
  // get the triggers IDs that will be used
  std::set<unsigned> triggers = GetTriggerIDs(FLAGS_triggers);
  
  // create centrality definition
  CentralityRun14 centrality;
  
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
    
    centrality.setEvent(header->GetRunId(), header->GetReferenceMultiplicity(),
                        header->GetZdcCoincidenceRate(), header->GetPrimaryVertexZ());
    int cent = centrality.centrality9();
    
    if (cent < 0)
      continue;
    
    double vz = header->GetPrimaryVertexZ();
    
    // get tracks & towers
    TList* tracks = reader->GetListOfSelectedTracks();
    TIter nextTrack(tracks);
    
    while(TStarJetPicoPrimaryTrack* track = (TStarJetPicoPrimaryTrack*) nextTrack()) {
      if (fabs(track->GetEta()) > 1.0 || track->GetPt() < FLAGS_ptmin)
        continue;
      
      if (vz < -10) {
        if (track->GetPhi() > 0) {
          pt_pos_near[cent]->Fill(track->GetPt());
          nhit_pos_near[cent]->Fill(track->GetNOfFittedHits());
          dca_pos_near[cent]->Fill(track->GetDCA());
        }
        else {
          pt_neg_near[cent]->Fill(track->GetPt());
          nhit_neg_near[cent]->Fill(track->GetNOfFittedHits());
          dca_neg_near[cent]->Fill(track->GetDCA());
        }
      }
      else if (vz > 10) {
        if (track->GetPhi() > 0) {
          pt_pos_far[cent]->Fill(track->GetPt());
          nhit_pos_far[cent]->Fill(track->GetNOfFittedHits());
          dca_pos_far[cent]->Fill(track->GetDCA());
        }
        else {
          pt_neg_far[cent]->Fill(track->GetPt());
          nhit_neg_far[cent]->Fill(track->GetNOfFittedHits());
          dca_neg_far[cent]->Fill(track->GetDCA());
        }
      }
    }
  }
  
    // print results
    // create our histogram and canvas options
    histogramOpts hopts;
    canvasOpts copts;
    canvasOpts coptslogz;
    coptslogz.log_z = true;
    canvasOpts coptslogy;
    coptslogy.log_y = true;
    canvasOpts cOptsBottomLeg;
    cOptsBottomLeg.leg_upper_bound = 0.4;
    cOptsBottomLeg.leg_lower_bound = 0.18;
    cOptsBottomLeg.leg_right_bound = 0.9;
    cOptsBottomLeg.leg_left_bound = 0.7;
    canvasOpts cOptsBottomLeftLeg;
    cOptsBottomLeftLeg.leg_upper_bound = 0.4;
    cOptsBottomLeftLeg.leg_lower_bound = 0.18;
    cOptsBottomLeftLeg.leg_right_bound = 0.18;
    cOptsBottomLeftLeg.leg_left_bound = 0.4;
    canvasOpts cOptsBottomLeftLegLogy;
    cOptsBottomLeftLegLogy.log_y = true;
    cOptsBottomLeftLegLogy.leg_upper_bound = 0.4;
    cOptsBottomLeftLegLogy.leg_lower_bound = 0.18;
    cOptsBottomLeftLegLogy.leg_right_bound = 0.18;
    cOptsBottomLeftLegLogy.leg_left_bound = 0.4;
    canvasOpts cOptsTopLeftLeg;
    cOptsTopLeftLeg.leg_right_bound = 0.18;
    cOptsTopLeftLeg.leg_left_bound = 0.4;
  
  std::vector<std::string> centrality_string{"0-5%", "5-10%", "10-20%",
    "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%"};
  for (int i = 0; i <= pt_neg_near.size(); ++i) {
    pt_neg_near[i]->Scale(1.0 / pt_neg_near[i]->Integral());
    pt_neg_far[i]->Scale(1.0 / pt_neg_far[i]->Integral());
    pt_pos_near[i]->Scale(1.0 / pt_pos_near[i]->Integral());
    pt_pos_far[i]->Scale(1.0 / pt_pos_far[i]->Integral());
    
    dca_neg_near[i]->Scale(1.0 / dca_neg_near[i]->Integral());
    dca_neg_far[i]->Scale(1.0 / dca_neg_far[i]->Integral());
    dca_pos_near[i]->Scale(1.0 / dca_pos_near[i]->Integral());
    dca_pos_far[i]->Scale(1.0 / dca_pos_far[i]->Integral());
    
    nhit_neg_near[i]->Scale(1.0 / nhit_neg_near[i]->Integral());
    nhit_neg_far[i]->Scale(1.0 / nhit_neg_far[i]->Integral());
    nhit_pos_near[i]->Scale(1.0 / nhit_pos_near[i]->Integral());
    nhit_pos_far[i]->Scale(1.0 / nhit_pos_far[i]->Integral());
    Overlay1D(pt_neg_near[i], pt_neg_far[i], "vz < -10", "vz > 10",
              hopts, coptslogy, FLAGS_outdir, MakeString("ptneg", i), "", "p_{T}",
              "fraction", "|#phi| < 0, " + centrality_string[i]);
    Overlay1D(pt_pos_near[i], pt_pos_far[i], "vz < -10", "vz > 10",
              hopts, coptslogy, FLAGS_outdir, MakeString("ptpos", i), "", "p_{T}",
              "fraction", "|#phi| > 0, " + centrality_string[i]);
    
    Overlay1D(dca_neg_near[i], dca_neg_far[i], "vz < -10", "vz > 10",
              hopts, copts, FLAGS_outdir, MakeString("dcaneg", i), "", "DCA [cm]",
              "fraction", "|#phi| < 0, " + centrality_string[i]);
    Overlay1D(dca_pos_near[i], dca_pos_far[i], "vz < -10", "vz > 10",
              hopts, copts, FLAGS_outdir, MakeString("dcapos", i), "", "DCA [cm]",
              "fraction", "|#phi| > 0, " + centrality_string[i]);
    
    Overlay1D(nhit_neg_near[i], nhit_neg_far[i], "vz < -10", "vz > 10",
              hopts, cOptsTopLeftLeg, FLAGS_outdir, MakeString("nhitneg", i), "", "fitted hits",
              "fraction", "|#phi| < 0, " + centrality_string[i]);
    
    Overlay1D(nhit_pos_near[i], nhit_pos_far[i], "vz < -10", "vz > 10",
              hopts, cOptsTopLeftLeg, FLAGS_outdir, MakeString("nhitpos", i), "", "fitted hits",
              "fraction", "|#phi| > 0, " + centrality_string[i]);
  }
  
  return 0;
}

//void Overlay1D(H* h1,
//               H* h2,
//               std::string h1_title,
//               std::string h2_title,
//               histogramOpts hopts,
//               canvasOpts copts,
//               std::string output_loc,
//               std::string output_name,
//               std::string canvas_title,
//               std::string x_axis_label,
//               std::string y_axis_label,
//               std::string legend_title = "")
