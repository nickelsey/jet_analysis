#include "jet_analysis/util/common.hh"
#include "jet_analysis/util/arg_helper.hh"
#include "jet_analysis/util/string_util.hh"
#include "jet_analysis/util/root_print_routines.hh"
#include "jet_analysis/centrality/centrality_run14.hh"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TF1.h"
#include "TStyle.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "TCanvas.h"

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

// use glog and gflags
#include "gflags/gflags.h"
#include "glog/stl_logging.h"
#include "glog/logging.h"


DEFINE_string(dataOld, "y7_aj_old.root", "one input root file");
DEFINE_string(dataNew, "y7_aj_new.root", "the other input root file");
DEFINE_string(tree, "LEAD_INIT_R_0.4_alg_2_pt_20_const_eta_1_const_pt_2_MATCH_R_0.4_alg_2_pt_0_const_eta_1_const_pt_0.2_SUB_INIT_R_0.4_alg_2_pt_10_const_eta_1_const_pt_2_MATCH_R_0.4_alg_2_pt_0_const_eta_1_const_pt_0.2", "input tree name");
DEFINE_string(outdir, "tmp", "output directory");

struct event {
  int runid;
  int eventid;
  int grefmult;
  int npart;
  TLorentzVector jl;
  TLorentzVector js;
  TLorentzVector jlm;
  TLorentzVector jsm;
  
};


typedef unordered_map<pair<int, int>, event, PairHash> event_map;

bool ReadInTree(string filename, string treename, event_map& map) {
  
  TFile input_file(filename.c_str(), "READ");
  if (input_file.IsOpen()) {
    TTreeReader reader(treename.c_str(), &input_file);
    TTreeReaderValue<int> eventid(reader, "eventid");
    TTreeReaderValue<int> runid(reader, "runid");
    TTreeReaderValue<int> gref(reader, "grefmult");
    TTreeReaderValue<int> nprim(reader, "npart");
    TTreeReaderValue<TLorentzVector> jl(reader, "jl");
    TTreeReaderValue<TLorentzVector> js(reader, "js");
    TTreeReaderValue<TLorentzVector> jlm(reader, "jlm");
    TTreeReaderValue<TLorentzVector> jsm(reader, "jsm");
    while (reader.Next()) {
      event evt;
      evt.runid = *runid;
      evt.eventid = *eventid;
      evt.grefmult = *gref;
      evt.npart = *nprim;
      evt.jl = *jl;
      evt.js = *js;
      evt.jlm = *jlm;
      evt.jsm = *jsm;
      
      map[{*runid, *eventid}] = evt;
    }
    
    return true;
  }
  else
  return false;
}

std::set<std::pair<int, int>> GetMatchedKeys(std::vector<event_map>& maps) {
  std::set<std::pair<int, int>> ret;
  
  for (auto& event : maps.front()) {
    
    auto& key = event.first;
    
    bool use_event = true;
    for (int i = 1; i < maps.size(); ++i) {
      if ((maps[i]).find(key) == (maps[i]).end()) {
        use_event = false;
        break;
      }
    }
    if (use_event) ret.insert(key);
  }
  return ret;
}

int main(int argc, char* argv[]) {
  
  string usage = "compares productions of run 14 on an event-by-event basis";
  
  gflags::SetUsageMessage(usage);
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);
  
  // set drawing preferences for histograms and graphs
  gStyle->SetOptStat(false);
  gStyle->SetOptFit(false);
  gStyle->SetLegendBorderSize(0);
  
  std::vector<event_map> evt_map(2);
  ReadInTree(FLAGS_dataNew, FLAGS_tree, evt_map[0]);
  ReadInTree(FLAGS_dataOld, FLAGS_tree, evt_map[1]);

  
  std::set<std::pair<int, int>> matched_keys = GetMatchedKeys(evt_map);
  
  // loop over new data
  for (auto& evt : evt_map[0]) {
    if (matched_keys.find(evt.first) == matched_keys.end()) {
      LOG(INFO) << "Did not find event: " << evt.first << " in old dataset";
      LOG(INFO) << " ---------------------------------------------------- ";
      LOG(INFO) << "grefmult: " << evt.second.grefmult;
      LOG(INFO) << "npart: " << evt.second.npart;
      LOG(INFO) << "lead pt: " << evt.second.jl.Pt();
      LOG(INFO) << "lead eta: " << evt.second.jl.Eta();
      LOG(INFO) << "lead phi: " << evt.second.jl.Phi();
      LOG(INFO) << "sublead pt: " << evt.second.js.Pt();
      LOG(INFO) << "sublead eta: " << evt.second.js.Eta();
      LOG(INFO) << "sublead phi: " << evt.second.js.Phi();
      double dphi = evt.second.jl.Phi() - evt.second.js.Phi();
      while (dphi < 0) {
        dphi += TMath::Pi() * 2.0;
      }
      while (dphi > TMath::Pi() * 2) {
        dphi -= TMath::Pi() * 2.0;
      }
      LOG(INFO) << "dphi: " << dphi;
      if (dphi < TMath::Pi() - 0.4)
        LOG(INFO) << "WTF";
    }
  }
  
  // loop over old data
  for (auto& evt : evt_map[1]) {
    if (matched_keys.find(evt.first) == matched_keys.end()) {
      LOG(INFO) << "Did not find event: " << evt.first << " in new dataset";
      LOG(INFO) << " ---------------------------------------------------- ";
      LOG(INFO) << "grefmult: " << evt.second.grefmult;
      LOG(INFO) << "npart: " << evt.second.npart;
      LOG(INFO) << "lead pt: " << evt.second.jl.Pt();
      LOG(INFO) << "lead eta: " << evt.second.jl.Eta();
      LOG(INFO) << "lead phi: " << evt.second.jl.Phi();
      LOG(INFO) << "sublead pt: " << evt.second.js.Pt();
      LOG(INFO) << "sublead eta: " << evt.second.js.Eta();
      LOG(INFO) << "sublead phi: " << evt.second.js.Phi();
      double dphi = evt.second.jl.Phi() - evt.second.js.Phi();
      while (dphi < 0) {
        dphi += TMath::Pi() * 2.0;
      }
      while (dphi > TMath::Pi() * 2) {
        dphi -= TMath::Pi() * 2.0;
      }
      LOG(INFO) << "dphi: " << dphi;
    }
  }
  
  for (auto& key : matched_keys) {
    auto& evt1 = evt_map[0][key];
    auto& evt2 = evt_map[1][key];
    
    if (evt1.npart != evt2.npart) {
      LOG(INFO) << "difference in nprim: " << evt1.npart << " " << evt2.npart;
    }
    if (fabs(evt1.jl.Pt() - evt2.jl.Pt()) > 0.0001) {
      LOG(INFO) << "pt difference???: " << evt1.jl.Pt() << " " << evt2.jl.Pt();
    }
  }
  
  return 0;
}
