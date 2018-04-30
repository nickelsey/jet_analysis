// print_results_with_systematics.cc

#include "jet_analysis/util/root_print_routines.hh"
#include "jet_analysis/util/arg_helper.hh"
#include "jet_analysis/util/string_util.hh"

#include "TROOT.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TKey.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TStyle.h"

#include <unordered_map>
#include <iostream>
#include <exception>

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

using std::string;

// adds to the map all TTrees that conform to
// the DijetWorker's naming convention
void GetTreesFromFile(const TFile& file, std::unordered_map<string, TTree*>& map) {
  TKey *key;
  TIter next(file.GetListOfKeys());
  
  while ((key = (TKey*) next())) {
    // check if its a TTree
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TTree"))
    continue;
    
    // check if its produced by the DijetWorker
    // (in a rather handwavy fashion)
    string tmp(key->GetName());
    if (tmp.find("LEAD_INIT") == string::npos ||
        tmp.find("SUB_INIT") == string::npos)
    continue;
    
    map.insert({tmp, (TTree*) key->ReadObj()});
  }
}

std::vector<TH1D*> SplitByCentrality(TH2D* h, std::vector<std::pair<int, int>>& centrality) {
  
  std::vector<TH1D*> ret;
  for (int i = 0; i < centrality.size(); ++i) {
    string name = string(h->GetName()) + std::to_string(i);
    TH1D* tmp = h->ProjectionY(name.c_str(),
                               centrality[i].first,
                               centrality[i].second);
    ret.push_back(tmp);
  }
  return ret;
}

std::vector<TH1D*> SplitByBin(TH2D* h) {
  
  std::vector<TH1D*> ret;
  for (int i = 1; i <= h->GetXaxis()->GetNbins(); ++i) {
    string name = string(h->GetName()) + std::to_string(i);
    TH1D* tmp = h->ProjectionY(name.c_str(), i, i);
    ret.push_back(tmp);
  }
  return ret;
}

std::vector<TH1D*> AddBins(std::vector<TH1D*> container, std::vector<std::pair<int, int>> bins) {
  
  std::vector<TH1D*> ret(bins.size(), nullptr);
  for (int i = 0; i < container.size(); ++i) {
    for (int j = 0; j < bins.size(); ++j) {
      if (i >=bins[j].first && i <= bins[j].second) {
        if (ret[i] == nullptr) {
          string name = string(h[i]->GetName()) + std::to_string(j);
          ret[i] = h->ProjectionY(name.c_str(), i, i);
        }
        else {
          ret[i]->Add(container[i]);
        }
      }
    }
  }
  return ret;
}

// command line option results
struct Options {
  string auau_file  = "";
  string pp_file    = "";
  string out_loc    = "results";
  string out_prefix = "";
  string sys_tow_m  = "";
  string sys_tow_p  = "";
  string sys_trk_m  = "";
  string sys_trk_p  = "";
  int cent_low      = 0;
  int cent_high     = 16;
};

int main(int argc, char* argv[]) {
  
  // set drawing preferences for histograms and graphs
  gStyle->SetOptStat(false);
  gStyle->SetOptFit(false);
  gStyle->SetOptTitle(1);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetHatchesSpacing(1.0);
  gStyle->SetHatchesLineWidth(2);
  
  // create histogram options
  histogramOpts hOpts;
  canvasOpts cOpts;
  
  Options opts;
  for (int i = 1; i < argc; ++i) {
    if (ParseStrFlag(string(argv[i]), "--auau", &opts.auau_file) ||
        ParseStrFlag(string(argv[i]), "--pp", &opts.pp_file) ||
        ParseStrFlag(string(argv[i]), "--outputDir", &opts.out_loc) ||
        ParseStrFlag(string(argv[i]), "--filePrefix", &opts.out_prefix) ||
        ParseStrFlag(string(argv[i]), "--towLow", &opts.sys_tow_m) ||
        ParseStrFlag(string(argv[i]), "--towHigh", &opts.sys_tow_p) ||
        ParseStrFlag(string(argv[i]), "--trackLow", &opts.sys_trk_m) ||
        ParseStrFlag(string(argv[i]), "--trackHigh", &opts.sys_tow_p) ||
        ParseIntFlag(string(argv[i]), "--centLow", &opts.cent_low) ||
        ParseIntFlag(string(argv[i]), "--centHigh", &opts.cent_high)) continue;
    std::cerr << "Unknown command line option: " << argv[i] << std::endl;
    return 1;
  }
  
  // check to make sure we have valid inputs
  std::vector<string> inputs{opts.auau_file, opts.pp_file, opts.sys_trk_p, opts.sys_trk_m,
                             opts.sys_tow_p, opts.sys_tow_m};
  for (auto file : inputs) {
    if (!boost::filesystem::exists(opts.pp_file)) {
      std::cout << "input file " << file;
      std::cout << "doesn't exist: exiting" << std::endl;
      return 1;
    }
  }
  
  // define centrality bins
  std::vector<std::pair<int, int>> centrality_5{{0, 0}, {1, 1}, {2,3}, {4, 7}, {8, 11}, {12, 15}};
  std::vector<std::string> centrality_5_string{"0-5%", "5-10%", "10-20%", "20-40%", "40-60%", "60-80%"};
  
  // read in the files
  TFile auau_file(opts.auau_file.c_str(), "READ");
  TFile pp_file(opts.pp_file.c_str(), "READ");
  TFile tow_p_file(opts.sys_tow_p.c_str(), "READ");
  TFile tow_m_file(opts.sys_tow_m.c_str(), "READ");
  TFile trk_p_file(opts.sys_trk_p.c_str(), "READ");
  TFile trk_m_file(opts.sys_trk_m.c_str(), "READ");
  
  // now we'll get the trees from the files, ignoring any objects
  // in the file that don't conform to the naming conventions from
  // the DijetWorker. There are also coincidence histograms to save
  std::vector<string> keys;
  std::unordered_map<std::string, TTree*> auau_trees;
  std::unordered_map<std::string, TTree*> pp_trees;
  std::unordered_map<std::string, TTree*> tow_p_trees;
  std::unordered_map<std::string, TTree*> tow_m_trees;
  std::unordered_map<std::string, TTree*> trk_p_trees;
  std::unordered_map<std::string, TTree*> trk_m_trees;
  
  
  GetTreesFromFile(auau_file, auau_trees);
  GetTreesFromFile(pp_file, pp_trees);
  GetTreesFromFile(tow_p_file, tow_p_trees);
  GetTreesFromFile(tow_m_file, tow_m_trees);
  GetTreesFromFile(trk_p_file, trk_p_trees);
  GetTreesFromFile(trk_m_file, trk_m_trees);
  
  // match keys
  for (auto entry : auau_trees) {
    if (pp_trees.find(entry.first) != pp_trees.end()) {
      keys.push_back(entry.first);
    }
  }
  
  std::unordered_map<std::string, TH1D*> auau_lead_count;
  std::unordered_map<std::string, TH1D*> auau_sub_count;
  std::unordered_map<std::string, TH1D*> pp_lead_count;
  std::unordered_map<std::string, TH1D*> pp_sub_count;
  
  for (auto key : keys) {
    auau_lead_count.insert({key, (TH1D*) auau_file.Get(MakeString(key, "_lead_count").c_str())});
    auau_sub_count.insert({key, (TH1D*) auau_file.Get(MakeString(key, "_sublead_count").c_str())});
    pp_lead_count.insert({key, (TH1D*) pp_file.Get(MakeString(key, "_lead_count").c_str())});
    pp_sub_count.insert({key, (TH1D*) pp_file.Get(MakeString(key, "_sublead_count").c_str())});
  }
  
  // save all histograms so we can do comparisons
  // between different keys if we want
  std::unordered_map<string, TH2D*> auau_hard_lead_pt;
  std::unordered_map<string, TH2D*> auau_hard_sub_pt;
  std::unordered_map<string, TH2D*> auau_match_lead_pt;
  std::unordered_map<string, TH2D*> auau_match_sub_pt;
  std::unordered_map<string, TH2D*> pp_hard_lead_pt;
  std::unordered_map<string, TH2D*> pp_hard_sub_pt;
  std::unordered_map<string, TH2D*> pp_match_lead_pt;
  std::unordered_map<string, TH2D*> pp_match_sub_pt;
  std::unordered_map<string, TH2D*> auau_hard_aj;
  std::unordered_map<string, TH2D*> auau_match_aj;
  std::unordered_map<string, TH2D*> pp_hard_aj;
  std::unordered_map<string, TH2D*> pp_match_aj;
  std::unordered_map<string, TH2D*> auau_dphi;
  std::unordered_map<string, TH2D*> pp_dphi;
  std::unordered_map<string, TH2D*> auau_hard_lead_rp;
  
  std::unordered_map<string, std::vector<TH1D*>> auau_hard_lead_pt_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_hard_sub_pt_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_match_lead_pt_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_match_sub_pt_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_hard_lead_pt_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_hard_sub_pt_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_match_lead_pt_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_match_sub_pt_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_hard_aj_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_match_aj_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_hard_aj_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_match_aj_cent;
  
  std::unordered_map<string, TH2D*> auau_hard_lead_rho;
  std::unordered_map<string, TH2D*> auau_hard_sub_rho;
  std::unordered_map<string, TH2D*> auau_match_lead_rho;
  std::unordered_map<string, TH2D*> auau_match_sub_rho;
  std::unordered_map<string, TH2D*> pp_hard_lead_rho;
  std::unordered_map<string, TH2D*> pp_hard_sub_rho;
  std::unordered_map<string, TH2D*> pp_match_lead_rho;
  std::unordered_map<string, TH2D*> pp_match_sub_rho;
  
  std::unordered_map<string, TH2D*> auau_hard_lead_sig;
  std::unordered_map<string, TH2D*> auau_hard_sub_sig;
  std::unordered_map<string, TH2D*> auau_match_lead_sig;
  std::unordered_map<string, TH2D*> auau_match_sub_sig;
  std::unordered_map<string, TH2D*> pp_hard_lead_sig;
  std::unordered_map<string, TH2D*> pp_hard_sub_sig;
  std::unordered_map<string, TH2D*> pp_match_lead_sig;
  std::unordered_map<string, TH2D*> pp_match_sub_sig;
  
  std::unordered_map<string, std::vector<TH1D*>> auau_hard_lead_rho_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_hard_sub_rho_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_match_lead_rho_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_match_sub_rho_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_hard_lead_rho_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_hard_sub_rho_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_match_lead_rho_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_match_sub_rho_cent;
  
  std::unordered_map<string, std::vector<TH1D*>> auau_hard_lead_sig_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_hard_sub_sig_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_match_lead_sig_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_match_sub_sig_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_hard_lead_sig_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_hard_sub_sig_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_match_lead_sig_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_match_sub_sig_cent;
  
  std::unordered_map<string, TH2D*> auau_npart;
  std::unordered_map<string, TH2D*> pp_npart;
  
  std::unordered_map<string, std::vector<TH1D*>> auau_npart_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_npart_cent;
  
  
  // count which key we are on
  int entry = -1;
  // loop over all matched trees
  for (auto key : keys) {
    
    //count which key we're on
    entry++;
    
    // make output directory
    string out_loc = opts.out_loc + "/" + key;
    boost::filesystem::path dir(out_loc.c_str());
    boost::filesystem::create_directories(dir);
    
    // key prefix for names so histograms don't get confused
    std::string key_prefix = "key_" + std::to_string(entry) + "_";
    
    // and create the file name prefix
    string file_prefix = out_loc + "/" + opts.out_prefix;
    
    // create coincidence measurement
    TH1D* auau_coincidence1 = auau_lead_count[key];
    TH1D* auau_coincidence2 = auau_sub_count[key];
    TH1D* pp_coincidence1 = pp_lead_count[key];
    TH1D* pp_coincidence2 = pp_sub_count[key];
    
    auau_coincidence1->RebinX(40);
    auau_coincidence2->RebinX(40);
    auau_coincidence2->Divide(auau_coincidence1);
    pp_coincidence1->RebinX(40);
    pp_coincidence2->RebinX(40);
    pp_coincidence2->Divide(pp_coincidence1);
    std::cout <<"pp coincidence: " << pp_coincidence2->GetBinContent(1) << std::endl;
    Overlay1D(auau_coincidence2, pp_coincidence2, "Au+Au", "P+P w/ efficiency",
              hOpts, cOpts, file_prefix, "dijet_coincidence", "refmult", "fraction", "", "Coincidence Rate");

    // process the trees
    
    TTree* auau_tree = auau_trees[key];
    TTree* pp_tree = pp_trees[key];
    
    TTreeReader auau_reader(auau_tree);
    TTreeReader pp_reader(pp_tree);
    
    // create readervalues for auau first
    TTreeReaderValue<int> auau_runid(auau_reader, "runid");
    TTreeReaderValue<int> auau_eventid(auau_reader, "eventid");
    TTreeReaderValue<double> auau_vz(auau_reader, "vz");
    TTreeReaderValue<int> auau_refmult(auau_reader, "refmult");
    TTreeReaderValue<int> auau_grefmult(auau_reader, "grefmult");
    TTreeReaderValue<double> auau_refmultcorr(auau_reader, "refmultcorr");
    TTreeReaderValue<double> auau_grefmultcorr(auau_reader, "grefmultcorr");
    TTreeReaderValue<int> auau_cent(auau_reader, "cent");
    TTreeReaderValue<double> auau_zdcrate(auau_reader, "zdcrate");
    TTreeReaderValue<double> auau_rp(auau_reader, "rp");
    TTreeReaderValue<int> auau_nglobal(auau_reader, "nglobal");
    TTreeReaderValue<int> auau_nprt(auau_reader, "npart");
    TTreeReaderValue<TLorentzVector> auau_jl(auau_reader, "jl");
    TTreeReaderValue<TLorentzVector> auau_js(auau_reader, "js");
    TTreeReaderValue<TLorentzVector> auau_jlm(auau_reader, "jlm");
    TTreeReaderValue<TLorentzVector> auau_jsm(auau_reader, "jsm");
    TTreeReaderValue<int> auau_jlconst(auau_reader, "jlconst");
    TTreeReaderValue<double> auau_jlrho(auau_reader, "jlrho");
    TTreeReaderValue<double> auau_jlsig(auau_reader, "jlsig");
    TTreeReaderValue<int> auau_jlmconst(auau_reader, "jlmconst");
    TTreeReaderValue<double> auau_jlmrho(auau_reader, "jlmrho");
    TTreeReaderValue<double> auau_jlmsig(auau_reader, "jlmsig");
    TTreeReaderValue<int> auau_jsconst(auau_reader, "jsconst");
    TTreeReaderValue<double> auau_jsrho(auau_reader, "jsrho");
    TTreeReaderValue<double> auau_jssig(auau_reader, "jssig");
    TTreeReaderValue<int> auau_jsmconst(auau_reader, "jsmconst");
    TTreeReaderValue<double> auau_jsmrho(auau_reader, "jsmrho");
    TTreeReaderValue<double> auau_jsmsig(auau_reader, "jsmsig");
    
    TTreeReaderValue<int> pp_runid(pp_reader, "runid");
    TTreeReaderValue<int> pp_eventid(pp_reader, "eventid");
    TTreeReaderValue<double> pp_vz(pp_reader, "vz");
    TTreeReaderValue<int> pp_refmult(pp_reader, "refmult");
    TTreeReaderValue<int> pp_grefmult(pp_reader, "grefmult");
    TTreeReaderValue<double> pp_refmultcorr(pp_reader, "refmultcorr");
    TTreeReaderValue<double> pp_grefmultcorr(pp_reader, "grefmultcorr");
    TTreeReaderValue<int> pp_cent(pp_reader, "cent");
    TTreeReaderValue<double> pp_zdcrate(pp_reader, "zdcrate");
    TTreeReaderValue<double> pp_rp(pp_reader, "rp");
    TTreeReaderValue<int> pp_nglobal(pp_reader, "nglobal");
    TTreeReaderValue<int> pp_nprt(pp_reader, "npart");
    TTreeReaderValue<TLorentzVector> pp_jl(pp_reader, "jl");
    TTreeReaderValue<TLorentzVector> pp_js(pp_reader, "js");
    TTreeReaderValue<TLorentzVector> pp_jlm(pp_reader, "jlm");
    TTreeReaderValue<TLorentzVector> pp_jsm(pp_reader, "jsm");
    TTreeReaderValue<int> pp_jlconst(pp_reader, "jlconst");
    TTreeReaderValue<double> pp_jlrho(pp_reader, "jlrho");
    TTreeReaderValue<double> pp_jlsig(pp_reader, "jlsig");
    TTreeReaderValue<int> pp_jlmconst(pp_reader, "jlmconst");
    TTreeReaderValue<double> pp_jlmrho(pp_reader, "jlmrho");
    TTreeReaderValue<double> pp_jlmsig(pp_reader, "jlmsig");
    TTreeReaderValue<int> pp_jsconst(pp_reader, "jsconst");
    TTreeReaderValue<double> pp_jsrho(pp_reader, "jsrho");
    TTreeReaderValue<double> pp_jssig(pp_reader, "jssig");
    TTreeReaderValue<int> pp_jsmconst(pp_reader, "jsmconst");
    TTreeReaderValue<double> pp_jsmrho(pp_reader, "jsmrho");
    TTreeReaderValue<double> pp_jsmsig(pp_reader, "jsmsig");
    
    // check if embedding was done for pp
    bool pp_embedded = false;
    if (pp_tree->GetBranch("embed_eventid"))
      pp_embedded = true;
    // build all embedding branches (only used if pp_embedded is true)
    TTreeReaderValue<int> embed_eventid(pp_reader, "embed_eventid");
    TTreeReaderValue<int> embed_runid(pp_reader, "embed_runid");
    TTreeReaderValue<int> embed_refmult(pp_reader, "embed_refmult");
    TTreeReaderValue<int> embed_grefmult(pp_reader, "embed_grefmult");
    TTreeReaderValue<int> embed_nprt(pp_reader, "embed_npart");
    TTreeReaderValue<double> embed_refmultcorr(pp_reader, "embed_refmultcorr");
    TTreeReaderValue<double> embed_grefmultcorr(pp_reader, "embed_grefmultcorr");
    TTreeReaderValue<int> embed_cent(pp_reader, "embed_cent");
    TTreeReaderValue<double> embed_rp(pp_reader, "embed_rp");
    TTreeReaderValue<double> embed_zdcrate(pp_reader, "embed_zdcrate");
    TTreeReaderValue<double> embed_vz(pp_reader, "embed_vz");
    
    // -----------------
    // define histograms
    // -----------------
    
    // jet pt
    TH2D* h_auau_hard_lead_pt = new TH2D(MakeString(key_prefix, "auauhardleadpt").c_str(),
                                         "p_{T}", 16, 0, 16, 100, 0, 100);
    TH2D* h_auau_hard_sub_pt = new TH2D(MakeString(key_prefix, "auauhardsubpt").c_str(),
                                        "p_{T}", 16, 0, 16, 100, 0, 100);
    TH2D* h_auau_match_lead_pt = new TH2D(MakeString(key_prefix, "auaumatchleadpt").c_str(),
                                          "p_{T}", 16, 0, 16, 100, 0, 100);
    TH2D* h_auau_match_sub_pt = new TH2D(MakeString(key_prefix, "auaumatchsubpt").c_str(),
                                         "p_{T}", 16, 0, 16, 100, 0, 100);
    TH2D* h_pp_hard_lead_pt = new TH2D(MakeString(key_prefix, "pphardleadpt").c_str(),
                                       "p_{T}", 16, 0, 16, 100, 0, 100);
    TH2D* h_pp_hard_sub_pt = new TH2D(MakeString(key_prefix, "pphardsubpt").c_str(),
                                      "p_{T}", 16, 0, 16, 100, 0, 100);
    TH2D* h_pp_match_lead_pt = new TH2D(MakeString(key_prefix, "ppmatchleadpt").c_str(),
                                        "p_{T}", 16, 0, 16, 100, 0, 100);
    TH2D* h_pp_match_sub_pt = new TH2D(MakeString(key_prefix, "ppmatchsubpt").c_str(),
                                       "p_{T}", 16, 0, 16, 100, 0, 100);
    
    // aj
    TH2D* h_auau_hard_aj = new TH2D(MakeString(key_prefix, "auauhardaj").c_str(),
                                    "A_{J}", 16, 0, 16, 30, 0, 0.9);
    TH2D* h_auau_match_aj = new TH2D(MakeString(key_prefix, "auaumatchaj").c_str(),
                                     "A_{J}", 16, 0, 16, 30, 0, 0.9);
    TH2D* h_pp_hard_aj = new TH2D(MakeString(key_prefix, "pphardaj").c_str(),
                                  "A_{J}", 16, 0, 16, 30, 0, 0.9);
    TH2D* h_pp_match_aj = new TH2D(MakeString(key_prefix, "ppmatchaj").c_str(),
                                   "A_{J}", 16, 0, 16, 30, 0, 0.9);
    
    // dphi
    TH2D* h_auau_dphi = new TH2D(MakeString(key_prefix, "auaudphi").c_str(),
                                 "d#phi", 16, 0, 16, 100, 0, 2*TMath::Pi());
    TH2D* h_pp_dphi = new TH2D(MakeString(key_prefix, "ppdphi").c_str(),
                               "d#phi", 16, 0, 16, 100, 0, 2*TMath::Pi());
    
    // lead jet - reaction plane
    TH2D* h_auau_hard_lead_rp = new TH2D(MakeString(key_prefix, "auauhardleaddphi").c_str(),
                                         "phi - rp", 9, -0.5, 8.5, 100, 0, 2*TMath::Pi());
    
    // rho and sigma
    TH2D* h_auau_hard_lead_rho = new TH2D(MakeString(key_prefix, "auauhardleadrho").c_str(),
                                          "auau hard lead rho", 16, 0, 16, 100, 0, 100);
    TH2D* h_auau_hard_lead_sig = new TH2D(MakeString(key_prefix, "auauhardleadsig").c_str(),
                                          "auau hard lead sig", 16, 0, 16, 100, 0, 20);
    TH2D* h_auau_hard_sub_rho = new TH2D(MakeString(key_prefix, "auauhardsubrho").c_str(),
                                         "auau hard sub rho", 16, 0, 16, 100, 0, 100);
    TH2D* h_auau_hard_sub_sig = new TH2D(MakeString(key_prefix, "auauhardsubsig").c_str(),
                                         "auau hard sub sig", 16, 0, 16, 100, 0, 20);
    TH2D* h_auau_match_lead_rho = new TH2D(MakeString(key_prefix, "auaumatchleadrho").c_str(),
                                           "auau match lead rho", 16, 0, 16, 100, 0, 100);
    TH2D* h_auau_match_lead_sig = new TH2D(MakeString(key_prefix, "auaumatchleadsig").c_str(),
                                           "auau match lead sig", 16, 0, 16, 100, 0, 20);
    TH2D* h_auau_match_sub_rho = new TH2D(MakeString(key_prefix, "auaumatchsubrho").c_str(),
                                          "auau match sub rho", 16, 0, 16, 100, 0, 100);
    TH2D* h_auau_match_sub_sig = new TH2D(MakeString(key_prefix, "auaumatchsubsig").c_str(),
                                          "auau match sub sig", 16, 0, 16, 100, 0, 20);
    
    TH2D* h_pp_hard_lead_rho = new TH2D(MakeString(key_prefix, "pphardleadrho").c_str(),
                                        "pp hard lead rho", 16, 0, 16, 100, 0, 100);
    TH2D* h_pp_hard_lead_sig = new TH2D(MakeString(key_prefix, "pphardleadsig").c_str(),
                                        "pp hard lead sig", 16, 0, 16, 100, 0, 20);
    TH2D* h_pp_hard_sub_rho = new TH2D(MakeString(key_prefix, "pphardsubrho").c_str(),
                                       "pp hard sub rho", 16, 0, 16, 100, 0, 100);
    TH2D* h_pp_hard_sub_sig = new TH2D(MakeString(key_prefix, "pphardsubsig").c_str(),
                                       "pp hard sub sig", 16, 0, 16, 100, 0, 20);
    TH2D* h_pp_match_lead_rho = new TH2D(MakeString(key_prefix, "ppmatchleadrho").c_str(),
                                         "pp match lead rho", 16, 0, 16, 100, 0, 100);
    TH2D* h_pp_match_lead_sig = new TH2D(MakeString(key_prefix, "ppmatchleadsig").c_str(),
                                         "pp match lead sig", 16, 0, 16, 100, 0, 20);
    TH2D* h_pp_match_sub_rho = new TH2D(MakeString(key_prefix, "ppmatchsubrho").c_str(),
                                        "pp match sub rho", 16, 0, 16, 100, 0, 100);
    TH2D* h_pp_match_sub_sig = new TH2D(MakeString(key_prefix, "ppmatchsubsig").c_str(),
                                        "pp match sub sig", 16, 0, 16, 100, 0, 20);
    
    TH2D* h_auau_npart = new TH2D(MakeString(key_prefix, "auaunpart").c_str(), ";refmult;nPart",
                                  16, 0, 16, 100, 0.5, 2500.5);
    TH2D* h_pp_npart = new TH2D(MakeString(key_prefix, "ppnpart").c_str(), ";refmult;nPart",
                                16, 0, 16, 100, 0.5, 2500.5);
    
    // insert into the dictionaries
    auau_hard_lead_pt.insert({key, h_auau_hard_lead_pt});
    auau_hard_sub_pt.insert({key, h_auau_hard_lead_pt});
    auau_match_lead_pt.insert({key, h_auau_match_lead_pt});
    auau_match_sub_pt.insert({key, h_auau_match_lead_pt});
    pp_hard_lead_pt.insert({key, h_pp_hard_lead_pt});
    pp_hard_sub_pt.insert({key, h_pp_hard_sub_pt});
    pp_match_lead_pt.insert({key, h_pp_match_lead_pt});
    pp_match_sub_pt.insert({key, h_pp_match_sub_pt});
    auau_hard_aj.insert({key, h_auau_hard_aj});
    auau_match_aj.insert({key, h_auau_match_aj});
    pp_hard_aj.insert({key, h_pp_hard_aj});
    pp_match_aj.insert({key, h_pp_match_aj});
    auau_dphi.insert({key, h_auau_dphi});
    pp_dphi.insert({key, h_pp_dphi});
    auau_hard_lead_rp.insert({key, h_auau_hard_lead_rp});
    
    auau_hard_lead_rho.insert({key, h_auau_hard_lead_rho});
    auau_hard_lead_sig.insert({key, h_auau_hard_lead_sig});
    auau_hard_sub_rho.insert({key, h_auau_hard_sub_rho});
    auau_hard_sub_sig.insert({key, h_auau_hard_sub_sig});
    auau_match_lead_rho.insert({key, h_auau_match_lead_rho});
    auau_match_lead_sig.insert({key, h_auau_match_lead_sig});
    auau_match_sub_rho.insert({key, h_auau_match_sub_rho});
    auau_match_sub_sig.insert({key, h_auau_match_sub_sig});
    pp_hard_lead_rho.insert({key, h_pp_hard_lead_rho});
    pp_hard_lead_sig.insert({key, h_pp_hard_lead_sig});
    pp_hard_sub_rho.insert({key, h_pp_hard_sub_rho});
    pp_hard_sub_sig.insert({key, h_pp_hard_sub_sig});
    pp_match_lead_rho.insert({key, h_pp_match_lead_rho});
    pp_match_lead_sig.insert({key, h_pp_match_lead_sig});
    pp_match_sub_rho.insert({key, h_pp_match_sub_rho});
    pp_match_sub_sig.insert({key, h_pp_match_sub_sig});
    
    auau_npart.insert({key, h_auau_npart});
    pp_npart.insert({key, h_pp_npart});
    
    // loop over the data & fill histograms
    while (auau_reader.Next()) {
      
      // auau jet pt
      h_auau_hard_lead_pt->Fill(*auau_cent, (*auau_jl).Pt());
      h_auau_hard_lead_rho->Fill(*auau_cent, *auau_jlrho);
      h_auau_hard_lead_sig->Fill(*auau_cent, *auau_jlsig);
      h_auau_hard_sub_pt->Fill(*auau_cent, (*auau_js).Pt());
      h_auau_hard_sub_rho->Fill(*auau_cent, *auau_jsrho);
      h_auau_hard_sub_sig->Fill(*auau_cent, *auau_jssig);
      h_auau_match_lead_pt->Fill(*auau_cent,(*auau_jlm).Pt());
      h_auau_match_lead_rho->Fill(*auau_cent, *auau_jlmrho);
      h_auau_match_lead_sig->Fill(*auau_cent, *auau_jlmsig);
      h_auau_match_sub_pt->Fill(*auau_cent, (*auau_jsm).Pt());
      h_auau_match_sub_rho->Fill(*auau_cent, *auau_jsmrho);
      h_auau_match_sub_sig->Fill(*auau_cent, *auau_jsmsig);
      
      // auau Aj
      h_auau_hard_aj->Fill(*auau_cent,
                           fabs((*auau_jl).Pt() - (*auau_js).Pt())/((*auau_jl).Pt() + (*auau_js).Pt()));
      h_auau_match_aj->Fill(*auau_cent,
                            fabs((*auau_jlm).Pt() - (*auau_jsm).Pt())/((*auau_jlm).Pt() + (*auau_jsm).Pt()));
      
      
      // auau dphi
      double dphi = (*auau_jl).Phi() - (*auau_js).Phi();
      double dphi_rp = (*auau_jl).Phi() - *auau_rp;
      
      //rotate dphi to be in [0,2*pi]
      while (dphi < 0)
        dphi += 2 * TMath::Pi();
      while (dphi > 2.0 * TMath::Pi())
        dphi -= 2.0 * TMath::Pi();
      while (dphi_rp < 0)
        dphi_rp += 2.0 * TMath::Pi();
      while (dphi_rp > 2.0 * TMath::Pi())
        dphi_rp -= 2.0 * TMath::Pi();
      
      h_auau_dphi->Fill(*auau_cent, dphi);
      h_auau_hard_lead_rp->Fill(*auau_cent, dphi_rp);
      
      // and refmult/npart
      h_auau_npart->Fill(*auau_cent, *auau_nprt);
      
    }
    
    while (pp_reader.Next()) {
      if (pp_embedded) {
        
        // pp single jet pt
        h_pp_hard_lead_pt->Fill(*pp_cent, (*pp_jl).Pt());
        h_pp_hard_lead_rho->Fill(*pp_cent, *pp_jlrho);
        h_pp_hard_lead_sig->Fill(*pp_cent, *pp_jlsig);
        h_pp_hard_sub_pt->Fill(*pp_cent, (*pp_js).Pt());
        h_pp_hard_sub_rho->Fill(*pp_cent, *pp_jsrho);
        h_pp_hard_sub_sig->Fill(*pp_cent, *pp_jssig);
        h_pp_match_lead_pt->Fill(*pp_cent, (*pp_jlm).Pt());
        h_pp_match_lead_rho->Fill(*pp_cent, *pp_jlmrho);
        h_pp_match_lead_sig->Fill(*pp_cent, *pp_jlmsig);
        h_pp_match_sub_pt->Fill(*pp_cent, (*pp_jsm).Pt());
        h_pp_match_sub_rho->Fill(*pp_cent, *pp_jsmrho);
        h_pp_match_sub_sig->Fill(*pp_cent, *pp_jsmsig);
        
        // pp Aj
        h_pp_hard_aj->Fill(*pp_cent,
                           fabs((*pp_jl).Pt() - (*pp_js).Pt())/((*pp_jl).Pt() + (*pp_js).Pt()));
        h_pp_match_aj->Fill(*pp_cent,
                            fabs((*pp_jlm).Pt() - (*pp_jsm).Pt())/((*pp_jlm).Pt() + (*pp_jsm).Pt()));
        
        // pp dphi
        double dphi = (*pp_jl).Phi() - (*pp_js).Phi();
        //rotate dphi to be in [0,2*pi]
        while (dphi < 0)
          dphi += 2 * TMath::Pi();
        while (dphi > 2.0 * TMath::Pi())
          dphi -= 2.0 * TMath::Pi();
        
        h_pp_dphi->Fill(*pp_cent, dphi);
        
        // npart
        h_pp_npart->Fill(*pp_cent, *pp_nprt + *embed_nprt);
      }
      else {
        
        // pp single jet pt
        h_pp_hard_lead_pt->Fill(*pp_cent, (*pp_jl).Pt());
        h_pp_hard_lead_rho->Fill(*pp_cent, *pp_jlrho);
        h_pp_hard_lead_sig->Fill(*pp_cent, *pp_jlsig);
        h_pp_hard_sub_pt->Fill(*pp_cent, (*pp_js).Pt());
        h_pp_hard_sub_rho->Fill(*pp_cent, *pp_jsrho);
        h_pp_hard_sub_sig->Fill(*pp_cent, *pp_jssig);
        h_pp_match_lead_pt->Fill(*pp_cent, (*pp_jlm).Pt());
        h_pp_match_lead_rho->Fill(*pp_cent, *pp_jlmrho);
        h_pp_match_lead_sig->Fill(*pp_cent, *pp_jlmsig);
        h_pp_match_sub_pt->Fill(*pp_cent, (*pp_jsm).Pt());
        h_pp_match_sub_rho->Fill(*pp_cent, *pp_jsmrho);
        h_pp_match_sub_sig->Fill(*pp_cent, *pp_jsmsig);
        
        // pp Aj
        h_pp_hard_aj->Fill(*pp_cent,
                           ((*pp_jl).Pt() - (*pp_js).Pt())/((*pp_jl).Pt() + (*pp_js).Pt()));
        h_pp_match_aj->Fill(*pp_cent,
                            ((*pp_jlm).Pt() - (*pp_jsm).Pt())/((*pp_jlm).Pt() + (*pp_jsm).Pt()));
        
        // pp dphi
        double dphi = (*pp_jl).Phi() - (*pp_js).Phi();
        //rotate dphi to be in [0,2*pi]
        while (dphi < 0)
        dphi += 2 * TMath::Pi();
        while (dphi > 2.0 * TMath::Pi())
        dphi -= 2.0 * TMath::Pi();
        
        h_pp_dphi->Fill(*pp_cent, dphi);
        
        // npart
        h_pp_npart->Fill(*pp_cent, *pp_nprt);
      }
    }
    
    // print dphi
    h_auau_dphi->Scale(1.0/h_auau_dphi->Integral());
    h_pp_dphi->Scale(1.0/h_pp_dphi->Integral());
    
    Overlay1D((TH1D*)h_auau_dphi->ProjectionY(), (TH1D*)h_pp_dphi->ProjectionY(),
              "Au+Au d#phi lead-sub", "P+P d#phi lead-sub", hOpts, cOpts,
              out_loc, "auau_pp_dphi", "", "d#phi", "fraction", "");
    Print2DSimple(h_auau_hard_lead_rp, hOpts, cOpts, out_loc, "auau_dphi_rp", "AuAu d#phi lead jet - rp",
                  "Centrality", "d#phi(lead-rp)");
    std::vector<TH1D*> h_rp_by_cent = SplitByCentrality(h_auau_hard_lead_rp, centrality_5);
    for (int i = 0; i < h_rp_by_cent.size(); ++i) {
      h_rp_by_cent[i]->Scale(1.0/h_rp_by_cent[i]->Integral());
      h_rp_by_cent[i]->RebinX(4);
    }
    Overlay1D(h_rp_by_cent, centrality_5_string, hOpts, cOpts, out_loc, "auau_cent_rp", "",
              "d#phi", "fraction", "Centrality");
    Overlay1D(h_rp_by_cent[0], h_rp_by_cent[1], "0-5%", "5-10%", hOpts, cOpts, out_loc, "auau_cent_rp_restricted", "",
              "d#phi", "fraction", "Centrality");
    
    // extract nPart in centrality bins
    std::vector<TH1D*> h_auau_npart_spectra = SplitByCentrality(h_auau_npart, centrality_5);
    std::vector<TH1D*> h_pp_npart_spectra = SplitByCentrality(h_pp_npart, centrality_5);
    
    for (int i = 0; i < h_auau_npart_spectra.size(); ++i) {
      h_auau_npart_spectra[i]->Scale(1.0/h_auau_npart_spectra[i]->Integral());
      h_pp_npart_spectra[i]->Scale(1.0/h_pp_npart_spectra[i]->Integral());
    }
    
    Overlay1D(h_auau_npart_spectra, centrality_5_string, hOpts, cOpts, out_loc, "auau_npart_spec",
              "", "npart", "fraction", "Centrality");
    Overlay1D(h_pp_npart_spectra, centrality_5_string, hOpts, cOpts, out_loc, "pp_npart_spec",
              "", "npart", "fraction", "Centrality");
    
    // extract pt spectra in centrality bins
    std::vector<TH1D*> h_auau_hard_lead_pt_spectra = SplitByCentrality(h_auau_hard_lead_pt,
                                                                    centrality_5);
    std::vector<TH1D*> h_auau_hard_sub_pt_spectra = SplitByCentrality(h_auau_hard_sub_pt,
                                                                   centrality_5);
    std::vector<TH1D*> h_auau_match_lead_pt_spectra = SplitByCentrality(h_auau_match_lead_pt,
                                                                     centrality_5);
    std::vector<TH1D*> h_auau_match_sub_pt_spectra = SplitByCentrality(h_auau_match_sub_pt,
                                                                    centrality_5);
    std::vector<TH1D*> h_pp_hard_lead_pt_spectra = SplitByCentrality(h_pp_hard_lead_pt,
                                                                  centrality_5);
    std::vector<TH1D*> h_pp_hard_sub_pt_spectra = SplitByCentrality(h_pp_hard_sub_pt,
                                                                 centrality_5);
    std::vector<TH1D*> h_pp_match_lead_pt_spectra = SplitByCentrality(h_pp_match_lead_pt,
                                                                   centrality_5);
    std::vector<TH1D*> h_pp_match_sub_pt_spectra = SplitByCentrality(h_pp_match_sub_pt,
                                                                  centrality_5);
    
    // extract rho & sigma in centrality bins
    std::vector<TH1D*> h_auau_hard_lead_rho_spectra = SplitByCentrality(h_auau_hard_lead_rho, centrality_5);
    std::vector<TH1D*> h_auau_hard_lead_sig_spectra = SplitByCentrality(h_auau_hard_lead_sig, centrality_5);
    std::vector<TH1D*> h_auau_hard_sub_rho_spectra = SplitByCentrality(h_auau_hard_sub_rho, centrality_5);
    std::vector<TH1D*> h_auau_hard_sub_sig_spectra = SplitByCentrality(h_auau_hard_sub_sig, centrality_5);
    std::vector<TH1D*> h_auau_match_lead_rho_spectra = SplitByCentrality(h_auau_match_lead_rho, centrality_5);
    std::vector<TH1D*> h_auau_match_lead_sig_spectra = SplitByCentrality(h_auau_match_lead_sig, centrality_5);
    std::vector<TH1D*> h_auau_match_sub_rho_spectra = SplitByCentrality(h_auau_match_sub_rho, centrality_5);
    std::vector<TH1D*> h_auau_match_sub_sig_spectra = SplitByCentrality(h_auau_match_sub_sig, centrality_5);
    std::vector<TH1D*> h_pp_hard_lead_rho_spectra = SplitByCentrality(h_pp_hard_lead_rho, centrality_5);
    std::vector<TH1D*> h_pp_hard_lead_sig_spectra = SplitByCentrality(h_pp_hard_lead_sig, centrality_5);
    std::vector<TH1D*> h_pp_hard_sub_rho_spectra = SplitByCentrality(h_pp_hard_sub_rho, centrality_5);
    std::vector<TH1D*> h_pp_hard_sub_sig_spectra = SplitByCentrality(h_pp_hard_sub_sig, centrality_5);
    std::vector<TH1D*> h_pp_match_lead_rho_spectra = SplitByCentrality(h_pp_match_lead_rho, centrality_5);
    std::vector<TH1D*> h_pp_match_lead_sig_spectra = SplitByCentrality(h_pp_match_lead_sig, centrality_5);
    std::vector<TH1D*> h_pp_match_sub_rho_spectra = SplitByCentrality(h_pp_match_sub_rho, centrality_5);
    std::vector<TH1D*> h_pp_match_sub_sig_spectra = SplitByCentrality(h_pp_match_sub_sig, centrality_5);
    
    
    // normalize pt spectra & rho & sig
    for (int i = 0; i < h_auau_hard_lead_pt_spectra.size(); ++i) {
      h_auau_hard_lead_pt_spectra[i]->Scale(1.0/h_auau_hard_lead_pt_spectra[i]->Integral());
      h_auau_hard_sub_pt_spectra[i]->Scale(1.0/h_auau_hard_sub_pt_spectra[i]->Integral());
      h_auau_match_lead_pt_spectra[i]->Scale(1.0/h_auau_match_lead_pt_spectra[i]->Integral());
      h_auau_match_sub_pt_spectra[i]->Scale(1.0/h_auau_match_sub_pt_spectra[i]->Integral());
      h_pp_hard_lead_pt_spectra[i]->Scale(1.0/h_pp_hard_lead_pt_spectra[i]->Integral());
      h_pp_hard_sub_pt_spectra[i]->Scale(1.0/h_pp_hard_sub_pt_spectra[i]->Integral());
      h_pp_match_lead_pt_spectra[i]->Scale(1.0/h_pp_match_lead_pt_spectra[i]->Integral());
      h_pp_match_sub_pt_spectra[i]->Scale(1.0/h_pp_match_sub_pt_spectra[i]->Integral());
      
      h_auau_hard_lead_rho_spectra[i]->Scale(1.0/h_auau_hard_lead_rho_spectra[i]->Integral());
      h_auau_hard_lead_sig_spectra[i]->Scale(1.0/h_auau_hard_lead_sig_spectra[i]->Integral());
      h_auau_hard_sub_rho_spectra[i]->Scale(1.0/h_auau_hard_sub_rho_spectra[i]->Integral());
      h_auau_hard_sub_sig_spectra[i]->Scale(1.0/h_auau_hard_sub_sig_spectra[i]->Integral());
      h_auau_match_lead_rho_spectra[i]->Scale(1.0/h_auau_match_lead_rho_spectra[i]->Integral());
      h_auau_match_lead_sig_spectra[i]->Scale(1.0/h_auau_match_lead_sig_spectra[i]->Integral());
      h_auau_match_sub_rho_spectra[i]->Scale(1.0/h_auau_match_sub_rho_spectra[i]->Integral());
      h_auau_match_sub_sig_spectra[i]->Scale(1.0/h_auau_match_sub_sig_spectra[i]->Integral());
      
      h_pp_hard_lead_rho_spectra[i]->Scale(1.0/h_pp_hard_lead_rho_spectra[i]->Integral());
      h_pp_hard_lead_sig_spectra[i]->Scale(1.0/h_pp_hard_lead_sig_spectra[i]->Integral());
      h_pp_hard_sub_rho_spectra[i]->Scale(1.0/h_pp_hard_sub_rho_spectra[i]->Integral());
      h_pp_hard_sub_sig_spectra[i]->Scale(1.0/h_pp_hard_sub_sig_spectra[i]->Integral());
      h_pp_match_lead_rho_spectra[i]->Scale(1.0/h_pp_match_lead_rho_spectra[i]->Integral());
      h_pp_match_lead_sig_spectra[i]->Scale(1.0/h_pp_match_lead_sig_spectra[i]->Integral());
      h_pp_match_sub_rho_spectra[i]->Scale(1.0/h_pp_match_sub_rho_spectra[i]->Integral());
      h_pp_match_sub_sig_spectra[i]->Scale(1.0/h_pp_match_sub_sig_spectra[i]->Integral());
    }
    
    // print pt spectra
    Overlay1D(h_auau_hard_lead_pt_spectra, centrality_5_string, hOpts, cOpts, out_loc,
              "auau_hard_lead_pt", "", "p_{T}", "fraction", "Centrality");
    Overlay1D(h_auau_hard_sub_pt_spectra, centrality_5_string, hOpts, cOpts, out_loc,
              "auau_hard_sub_pt", "", "p_{T}", "fraction", "Centrality");
    Overlay1D(h_auau_match_lead_pt_spectra, centrality_5_string, hOpts, cOpts, out_loc,
              "auau_match_lead_pt", "", "p_{T}", "fraction", "Centrality");
    Overlay1D(h_auau_match_sub_pt_spectra, centrality_5_string, hOpts, cOpts, out_loc,
              "auau_match_sub_pt", "", "p_{T}", "fraction", "Centrality");
    Overlay1D(h_pp_hard_lead_pt_spectra, centrality_5_string, hOpts, cOpts, out_loc,
              "pp_hard_lead_pt", "", "p_{T}", "fraction", "Centrality");
    Overlay1D(h_pp_hard_sub_pt_spectra, centrality_5_string, hOpts, cOpts, out_loc,
              "pp_hard_sub_pt", "", "p_{T}", "fraction", "Centrality");
    Overlay1D(h_pp_match_lead_pt_spectra, centrality_5_string, hOpts, cOpts, out_loc,
              "pp_match_lead_pt", "", "p_{T}", "fraction", "Centrality");
    Overlay1D(h_pp_match_sub_pt_spectra, centrality_5_string, hOpts, cOpts, out_loc,
              "pp_match_sub_pt", "", "p_{T}", "fraction", "Centrality");
    
    // print rho & sig
    Overlay1D(h_auau_hard_lead_rho_spectra, centrality_5_string, hOpts, cOpts, out_loc,
              "auau_hard_lead_rho", "", "#rho", "fraction", "Centrality");
    Overlay1D(h_auau_hard_lead_sig_spectra, centrality_5_string, hOpts, cOpts, out_loc,
              "auau_hard_lead_sig", "", "#sigma", "fraction", "Centrality");
    Overlay1D(h_auau_hard_sub_rho_spectra, centrality_5_string, hOpts, cOpts, out_loc,
              "auau_hard_sub_rho", "", "#rho", "fraction", "Centrality");
    Overlay1D(h_auau_hard_sub_sig_spectra, centrality_5_string, hOpts, cOpts, out_loc,
              "auau_hard_sub_sig", "", "#sigma", "fraction", "Centrality");
    Overlay1D(h_auau_match_lead_rho_spectra, centrality_5_string, hOpts, cOpts, out_loc,
              "auau_match_lead_rho", "", "#rho", "fraction", "Centrality");
    Overlay1D(h_auau_match_lead_sig_spectra, centrality_5_string, hOpts, cOpts, out_loc,
              "auau_match_lead_sig", "", "#sigma", "fraction", "Centrality");
    Overlay1D(h_auau_match_sub_rho_spectra, centrality_5_string, hOpts, cOpts, out_loc,
              "auau_match_sub_rho", "", "#rho", "fraction", "Centrality");
    Overlay1D(h_auau_match_sub_sig_spectra, centrality_5_string, hOpts, cOpts, out_loc,
              "auau_match_sub_sig", "", "#sigma", "fraction", "Centrality");
    
    Overlay1D(h_pp_hard_lead_rho_spectra, centrality_5_string, hOpts, cOpts, out_loc,
              "pp_hard_lead_rho", "", "#rho", "fraction", "Centrality");
    Overlay1D(h_pp_hard_lead_sig_spectra, centrality_5_string, hOpts, cOpts, out_loc,
              "pp_hard_lead_sig", "", "#sigma", "fraction","Centrality");
    Overlay1D(h_pp_hard_sub_rho_spectra, centrality_5_string, hOpts, cOpts, out_loc,
              "pp_hard_sub_rho", "", "#rho", "fraction", "Centrality");
    Overlay1D(h_pp_hard_sub_sig_spectra, centrality_5_string, hOpts, cOpts, out_loc,
              "pp_hard_sub_sig", "", "#sigma", "fraction", "Centrality");
    Overlay1D(h_pp_match_lead_rho_spectra, centrality_5_string, hOpts, cOpts, out_loc,
              "pp_match_lead_rho", "", "#rho", "fraction", "Centrality");
    Overlay1D(h_pp_match_lead_sig_spectra, centrality_5_string, hOpts, cOpts, out_loc,
              "pp_match_lead_sig", "", "#sigma", "fraction", "Centrality");
    Overlay1D(h_pp_match_sub_rho_spectra, centrality_5_string, hOpts, cOpts, out_loc,
              "pp_match_sub_rho", "", "#rho", "fraction", "Centrality");
    Overlay1D(h_pp_match_sub_sig_spectra, centrality_5_string, hOpts, cOpts, out_loc,
              "pp_match_sub_sig", "", "#sigma", "fraction", "Centrality");
    
    // for Aj
    std::vector<TH1D*> h_auau_hard_aj_spectra = SplitByBin(h_auau_hard_aj);
    std::vector<TH1D*> h_auau_match_aj_spectra = SplitByBin(h_auau_match_aj);

    std::vector<TH1D*> h_pp_hard_aj_spectra = SplitByBin(h_pp_hard_aj);
    std::vector<TH1D*> h_pp_match_aj_spectra = SplitByBin(h_pp_match_aj);
    
    // we will weight pp by the relative fraction of events in auau bins
                                 
                                 
    for (int i = 0; i < h_auau_hard_aj_spectra.size(); ++i) {
      h_auau_hard_aj_spectra[i]->RebinX(2);
      h_auau_match_aj_spectra[i]->RebinX(2);
      h_pp_hard_aj_spectra[i]->RebinX(2);
      h_pp_match_aj_spectra[i]->RebinX(2);
      h_auau_hard_aj_spectra[i]->Scale(1.0/h_auau_hard_aj_spectra[i]->Integral());
      h_auau_match_aj_spectra[i]->Scale(1.0/h_auau_match_aj_spectra[i]->Integral());
      h_pp_hard_aj_spectra[i]->Scale(1.0/h_pp_hard_aj_spectra[i]->Integral());
      h_pp_match_aj_spectra[i]->Scale(1.0/h_pp_match_aj_spectra[i]->Integral());
    }
//
//    Overlay1D(h_auau_hard_aj_spectra, centrality_5_string, hOpts, cOpts, out_loc, "auau_hard_aj", "",
//              "A_{J}", "fraction", "Centrality");
//    Overlay1D(h_auau_match_aj_spectra, centrality_5_string, hOpts, cOpts, out_loc, "auau_match_aj", "",
//              "A_{J}", "fraction", "Centrality");
//    Overlay1D(h_pp_hard_aj_spectra, centrality_5_string, hOpts, cOpts, out_loc, "pp_hard_aj", "",
//              "A_{J}", "fraction", "Centrality");
//    Overlay1D(h_pp_match_aj_spectra, centrality_5_string, hOpts, cOpts, out_loc, "pp_match_aj", "",
//              "A_{J}", "fraction", "Centrality");
//
//    // add the containers to the dictionaries
//    auau_hard_lead_pt_cent.insert({key, h_auau_hard_lead_pt_spectra});
//    auau_hard_sub_pt_cent.insert({key, h_auau_hard_sub_pt_spectra});
//    auau_match_lead_pt_cent.insert({key, h_auau_match_lead_pt_spectra});
//    auau_match_sub_pt_cent.insert({key, h_auau_match_sub_pt_spectra});
//    pp_hard_lead_pt_cent.insert({key, h_pp_hard_lead_pt_spectra});
//    pp_hard_sub_pt_cent.insert({key, h_pp_hard_sub_pt_spectra});
//    pp_match_lead_pt_cent.insert({key, h_pp_match_lead_pt_spectra});
//    pp_match_sub_pt_cent.insert({key, h_pp_match_sub_pt_spectra});
//    auau_hard_aj_cent.insert({key, h_auau_hard_aj_spectra});
//    auau_match_aj_cent.insert({key, h_auau_match_aj_spectra});
//    auau_hard_aj_cent.insert({key, h_pp_hard_aj_spectra});
//    auau_match_aj_cent.insert({key, h_pp_match_aj_spectra});
//
//    auau_hard_lead_rho_cent.insert({key, h_auau_hard_lead_rho_spectra});
//    auau_hard_lead_sig_cent.insert({key, h_auau_hard_lead_sig_spectra});
//    auau_hard_sub_rho_cent.insert({key, h_auau_hard_sub_rho_spectra});
//    auau_hard_sub_sig_cent.insert({key, h_auau_hard_sub_sig_spectra});
//    auau_match_lead_rho_cent.insert({key, h_auau_match_lead_rho_spectra});
//    auau_match_lead_sig_cent.insert({key, h_auau_match_lead_sig_spectra});
//    auau_match_sub_rho_cent.insert({key, h_auau_match_sub_rho_spectra});
//    auau_match_sub_sig_cent.insert({key, h_auau_match_sub_sig_spectra});
//
//    pp_hard_lead_rho_cent.insert({key, h_pp_hard_lead_rho_spectra});
//    pp_hard_lead_sig_cent.insert({key, h_pp_hard_lead_sig_spectra});
//    pp_hard_sub_rho_cent.insert({key, h_pp_hard_sub_rho_spectra});
//    pp_hard_sub_sig_cent.insert({key, h_pp_hard_sub_sig_spectra});
//    pp_match_lead_rho_cent.insert({key, h_pp_match_lead_rho_spectra});
//    pp_match_lead_sig_cent.insert({key, h_pp_match_lead_sig_spectra});
//    pp_match_sub_rho_cent.insert({key, h_pp_match_sub_rho_spectra});
//    pp_match_sub_sig_cent.insert({key, h_pp_match_sub_sig_spectra});
//
//    // now, printing some comparisons
//    // make a folder for each centrality
//    // make output directory
//    for (int i = 0; i < h_auau_hard_aj_spectra.size(); ++i) {
//      std::string out_loc_cent = out_loc + "/cent_" + std::to_string(i);
//      boost::filesystem::path dir(out_loc_cent.c_str());
//      boost::filesystem::create_directories(dir);
//
//      Overlay1D(h_auau_hard_aj_spectra[i], h_pp_hard_aj_spectra[i], "AuAu hard A_{J}", "PP hard A_{J}",
//                 hOpts, cOpts, out_loc_cent, "aj_hard", "", "A_{J}", "fraction", "A_{J}");
//      Overlay1D(h_auau_match_aj_spectra[i], h_pp_match_aj_spectra[i], "AuAu matched A_{J}", "PP matched A_{J}",
//                 hOpts, cOpts, out_loc_cent, "aj_match", "", "A_{J}", "fraction", "A_{J}");
//    }
  }
  
  return 0;
}

