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
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLatex.h"

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
                               centrality[i].first + 1,
                               centrality[i].second + 1);
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
        if (ret[j] == nullptr) {
          string name = string(container[i]->GetName()) + std::to_string(j);
          ret[j] = (TH1D*) container[i]->Clone(name.c_str());
        }
        else {
          ret[j]->Add(container[i]);
        }
      }
    }
  }
  return ret;
}

TGraphErrors* GetSystematic(TH1D* nom, TH1D* var1_a, TH1D* var1_b, TH1D* var2_a, TH1D* var2_b) {
  int nBins = nom->GetNbinsX();
  double x_[nBins];
  double y_[nBins];
  double x_err_[nBins];
  double y_err_[nBins];
  
  for (int i = 0; i < nBins; ++i) {
    x_[i] = nom->GetBinCenter(i+1);
    y_[i] = nom->GetBinContent(i+1);
    x_err_[i] = nom->GetXaxis()->GetBinWidth(1) / 2.0;
    double diff_var_1_a = fabs(nom->GetBinContent(i+1)  - var1_a->GetBinContent(i+1));
    double diff_var_1_b = fabs(nom->GetBinContent(i+1)  - var1_b->GetBinContent(i+1));
    double diff_var_2_a = fabs(nom->GetBinContent(i+1)  - var2_a->GetBinContent(i+1));
    double diff_var_2_b = fabs(nom->GetBinContent(i+1)  - var2_b->GetBinContent(i+1));
    double max_var_1 = (diff_var_1_a > diff_var_1_b ? diff_var_1_a : diff_var_1_a);
    double max_var_2 = (diff_var_2_a > diff_var_2_b ? diff_var_2_a : diff_var_2_b);
    y_err_[i] = sqrt(max_var_1 * max_var_1 + max_var_2 * max_var_2);
  }
  TGraphErrors* ret = new TGraphErrors(nBins, x_, y_, x_err_, y_err_);
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
        ParseStrFlag(string(argv[i]), "--trackHigh", &opts.sys_trk_p) ||
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
  //std::vector<std::pair<int, int>> centrality_5{{0, 0}, {1, 1}, {2,3}, {4, 7}, {8, 11}, {12, 15}};
  //std::vector<std::string> centrality_5_string{"0-5%", "5-10%", "10-20%", "20-40%", "40-60%", "60-80%"};
  
  std::vector<std::pair<int, int>> centrality_5{{0, 3}, {4, 7}, {10,13}};
  std::vector<std::string> centrality_5_string{"0-20%", "20-40%", "50-70%"};
  
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
  
  // for systematics
  std::unordered_map<string, TH2D*> pp_hard_aj_tow_p;
  std::unordered_map<string, TH2D*> pp_match_aj_tow_p;
  std::unordered_map<string, TH2D*> pp_hard_aj_tow_m;
  std::unordered_map<string, TH2D*> pp_match_aj_tow_m;
  std::unordered_map<string, TH2D*> pp_hard_aj_trk_p;
  std::unordered_map<string, TH2D*> pp_match_aj_trk_p;
  std::unordered_map<string, TH2D*> pp_hard_aj_trk_m;
  std::unordered_map<string, TH2D*> pp_match_aj_trk_m;
  
  std::unordered_map<string, std::vector<TH1D*>> pp_hard_aj_tow_p_spectra;
  std::unordered_map<string, std::vector<TH1D*>> pp_match_aj_tow_p_spectra;
  std::unordered_map<string, std::vector<TH1D*>> pp_hard_aj_tow_m_spectra;
  std::unordered_map<string, std::vector<TH1D*>> pp_match_aj_tow_m_spectra;
  std::unordered_map<string, std::vector<TH1D*>> pp_hard_aj_trk_p_spectra;
  std::unordered_map<string, std::vector<TH1D*>> pp_match_aj_trk_p_spectra;
  std::unordered_map<string, std::vector<TH1D*>> pp_hard_aj_trk_m_spectra;
  std::unordered_map<string, std::vector<TH1D*>> pp_match_aj_trk_m_spectra;
  
  
  
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
    
    TTree* tow_p_tree = tow_p_trees[key];
    TTree* tow_m_tree = tow_m_trees[key];
    TTree* trk_p_tree = trk_p_trees[key];
    TTree* trk_m_tree = trk_m_trees[key];
    
    TTreeReader auau_reader(auau_tree);
    TTreeReader pp_reader(pp_tree);
    
    TTreeReader tow_p_reader(tow_p_tree);
    TTreeReader tow_m_reader(tow_m_tree);
    TTreeReader trk_p_reader(trk_p_tree);
    TTreeReader trk_m_reader(trk_m_tree);
    
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
    
    // and systematics
    TTreeReaderValue<TLorentzVector> tow_p_jl(tow_p_reader, "jl");
    TTreeReaderValue<TLorentzVector> tow_p_js(tow_p_reader, "js");
    TTreeReaderValue<TLorentzVector> tow_p_jlm(tow_p_reader, "jlm");
    TTreeReaderValue<TLorentzVector> tow_p_jsm(tow_p_reader, "jsm");
    TTreeReaderValue<int> tow_p_cent(tow_p_reader, "cent");
    
    TTreeReaderValue<TLorentzVector> tow_m_jl(tow_m_reader, "jl");
    TTreeReaderValue<TLorentzVector> tow_m_js(tow_m_reader, "js");
    TTreeReaderValue<TLorentzVector> tow_m_jlm(tow_m_reader, "jlm");
    TTreeReaderValue<TLorentzVector> tow_m_jsm(tow_m_reader, "jsm");
    TTreeReaderValue<int> tow_m_cent(tow_m_reader, "cent");
    
    TTreeReaderValue<TLorentzVector> trk_p_jl(trk_p_reader, "jl");
    TTreeReaderValue<TLorentzVector> trk_p_js(trk_p_reader, "js");
    TTreeReaderValue<TLorentzVector> trk_p_jlm(trk_p_reader, "jlm");
    TTreeReaderValue<TLorentzVector> trk_p_jsm(trk_p_reader, "jsm");
    TTreeReaderValue<int> trk_p_cent(trk_p_reader, "cent");
    
    TTreeReaderValue<TLorentzVector> trk_m_jl(trk_m_reader, "jl");
    TTreeReaderValue<TLorentzVector> trk_m_js(trk_m_reader, "js");
    TTreeReaderValue<TLorentzVector> trk_m_jlm(trk_m_reader, "jlm");
    TTreeReaderValue<TLorentzVector> trk_m_jsm(trk_m_reader, "jsm");
    TTreeReaderValue<int> trk_m_cent(trk_m_reader, "cent");
    
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
    
    TH2D* h_tow_p_hard_aj = new TH2D(MakeString(key_prefix, "towphardaj").c_str(),
                                  "A_{J}", 16, 0, 16, 30, 0, 0.9);
    TH2D* h_tow_p_match_aj = new TH2D(MakeString(key_prefix, "towpmatchaj").c_str(),
                                   "A_{J}", 16, 0, 16, 30, 0, 0.9);
    TH2D* h_tow_m_hard_aj = new TH2D(MakeString(key_prefix, "towmhardaj").c_str(),
                                     "A_{J}", 16, 0, 16, 30, 0, 0.9);
    TH2D* h_tow_m_match_aj = new TH2D(MakeString(key_prefix, "towmmatchaj").c_str(),
                                      "A_{J}", 16, 0, 16, 30, 0, 0.9);
    TH2D* h_trk_p_hard_aj = new TH2D(MakeString(key_prefix, "trkphardaj").c_str(),
                                     "A_{J}", 16, 0, 16, 30, 0, 0.9);
    TH2D* h_trk_p_match_aj = new TH2D(MakeString(key_prefix, "trkpmatchaj").c_str(),
                                      "A_{J}", 16, 0, 16, 30, 0, 0.9);
    TH2D* h_trk_m_hard_aj = new TH2D(MakeString(key_prefix, "trkmhardaj").c_str(),
                                     "A_{J}", 16, 0, 16, 30, 0, 0.9);
    TH2D* h_trk_m_match_aj = new TH2D(MakeString(key_prefix, "trkmmatchaj").c_str(),
                                      "A_{J}", 16, 0, 16, 30, 0, 0.9);
    
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
    
    pp_hard_aj_tow_p.insert({key, h_tow_p_hard_aj});
    pp_match_aj_tow_p.insert({key, h_tow_p_match_aj});
    pp_hard_aj_tow_m.insert({key, h_tow_m_hard_aj});
    pp_match_aj_tow_m.insert({key, h_tow_m_match_aj});
    pp_hard_aj_trk_p.insert({key, h_trk_p_hard_aj});
    pp_match_aj_trk_p.insert({key, h_trk_p_match_aj});
    pp_hard_aj_trk_m.insert({key, h_trk_m_hard_aj});
    pp_match_aj_trk_m.insert({key, h_trk_m_match_aj});
    
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
        h_pp_npart->Fill(*pp_cent, *pp_nprt);
      }
    }
    
    // next fill the systematics
    while (tow_p_reader.Next()) {
      h_tow_p_hard_aj->Fill(*tow_p_cent,
                            fabs((*tow_p_jl).Pt() - (*tow_p_js).Pt()) / ((*tow_p_jl).Pt() + (*tow_p_js).Pt()));
      h_tow_p_match_aj->Fill(*tow_p_cent,
                            fabs((*tow_p_jlm).Pt() - (*tow_p_jsm).Pt()) / ((*tow_p_jlm).Pt() + (*tow_p_jsm).Pt()));
    }
    
    while (tow_m_reader.Next()) {
      h_tow_m_hard_aj->Fill(*tow_m_cent,
                            fabs((*tow_m_jl).Pt() - (*tow_m_js).Pt()) / ((*tow_m_jl).Pt() + (*tow_m_js).Pt()));
      h_tow_m_match_aj->Fill(*tow_m_cent,
                            fabs((*tow_m_jlm).Pt() - (*tow_m_jsm).Pt()) / ((*tow_m_jlm).Pt() + (*tow_m_jsm).Pt()));
    }
    
    while (trk_p_reader.Next()) {
      h_trk_p_hard_aj->Fill(*trk_p_cent,
                            fabs((*trk_p_jl).Pt() - (*trk_p_js).Pt()) / ((*trk_p_jl).Pt() + (*trk_p_js).Pt()));
      h_trk_p_match_aj->Fill(*trk_p_cent,
                            fabs((*trk_p_jlm).Pt() - (*trk_p_jsm).Pt()) / ((*trk_p_jlm).Pt() + (*trk_p_jsm).Pt()));
    }
    
    while (trk_m_reader.Next()) {
      h_trk_m_hard_aj->Fill(*trk_m_cent,
                            fabs((*trk_m_jl).Pt() - (*trk_m_js).Pt()) / ((*trk_m_jl).Pt() + (*trk_m_js).Pt()));
      h_trk_m_match_aj->Fill(*trk_m_cent,
                            fabs((*trk_m_jlm).Pt() - (*trk_m_jsm).Pt()) / ((*trk_m_jlm).Pt() + (*trk_m_jsm).Pt()));
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
    std::vector<TH1D*> h_auau_hard_aj_spectra_bin = SplitByBin(h_auau_hard_aj);
    std::vector<TH1D*> h_auau_match_aj_spectra_bin = SplitByBin(h_auau_match_aj);

    std::vector<TH1D*> h_pp_hard_aj_spectra_bin = SplitByBin(h_pp_hard_aj);
    std::vector<TH1D*> h_pp_match_aj_spectra_bin = SplitByBin(h_pp_match_aj);
    
    std::vector<TH1D*> h_tow_p_hard_aj_spectra_bin = SplitByBin(h_tow_p_hard_aj);
    std::vector<TH1D*> h_tow_p_match_aj_spectra_bin = SplitByBin(h_tow_p_match_aj);
    std::vector<TH1D*> h_tow_m_hard_aj_spectra_bin = SplitByBin(h_tow_m_hard_aj);
    std::vector<TH1D*> h_tow_m_match_aj_spectra_bin = SplitByBin(h_tow_m_match_aj);
    
    std::vector<TH1D*> h_trk_p_hard_aj_spectra_bin = SplitByBin(h_trk_p_hard_aj);
    std::vector<TH1D*> h_trk_p_match_aj_spectra_bin = SplitByBin(h_trk_p_match_aj);
    std::vector<TH1D*> h_trk_m_hard_aj_spectra_bin = SplitByBin(h_trk_m_hard_aj);
    std::vector<TH1D*> h_trk_m_match_aj_spectra_bin = SplitByBin(h_trk_m_match_aj);
    
    // we will weight pp by the relative fraction of events in auau bins
    std::vector<double> auau_weights;
    double sum;
    for (auto h : h_auau_hard_aj_spectra_bin) {
      auau_weights.push_back(h->Integral());
      sum += h->Integral();
    }
    
    for (auto& entry : auau_weights) {
      entry /= sum;
    }
    
    for (int i = 0; i < h_auau_hard_aj_spectra_bin.size(); ++i) {
      
      h_auau_hard_aj_spectra_bin[i]->RebinX(2);
      h_auau_match_aj_spectra_bin[i]->RebinX(2);
      h_pp_hard_aj_spectra_bin[i]->RebinX(2);
      h_pp_match_aj_spectra_bin[i]->RebinX(2);
      
      h_auau_hard_aj_spectra_bin[i]->Scale(1.0/h_auau_hard_aj_spectra_bin[i]->Integral());
      h_auau_match_aj_spectra_bin[i]->Scale(1.0/h_auau_match_aj_spectra_bin[i]->Integral());
      h_pp_hard_aj_spectra_bin[i]->Scale(auau_weights[i]/h_pp_hard_aj_spectra_bin[i]->Integral());
      h_pp_match_aj_spectra_bin[i]->Scale(auau_weights[i]/h_pp_match_aj_spectra_bin[i]->Integral());
      
      h_tow_p_hard_aj_spectra_bin[i]->RebinX(2);
      h_tow_p_match_aj_spectra_bin[i]->RebinX(2);
      h_tow_m_hard_aj_spectra_bin[i]->RebinX(2);
      h_tow_m_match_aj_spectra_bin[i]->RebinX(2);
      
      h_trk_p_hard_aj_spectra_bin[i]->RebinX(2);
      h_trk_p_match_aj_spectra_bin[i]->RebinX(2);
      h_trk_m_hard_aj_spectra_bin[i]->RebinX(2);
      h_trk_m_match_aj_spectra_bin[i]->RebinX(2);
      
      h_tow_p_hard_aj_spectra_bin[i]->Scale(auau_weights[i]/h_tow_p_hard_aj_spectra_bin[i]->Integral());
      h_tow_p_match_aj_spectra_bin[i]->Scale(auau_weights[i]/h_tow_p_match_aj_spectra_bin[i]->Integral());
      h_tow_m_hard_aj_spectra_bin[i]->Scale(auau_weights[i]/h_tow_m_hard_aj_spectra_bin[i]->Integral());
      h_tow_m_match_aj_spectra_bin[i]->Scale(auau_weights[i]/h_tow_m_match_aj_spectra_bin[i]->Integral());
      
      h_trk_p_hard_aj_spectra_bin[i]->Scale(auau_weights[i]/h_trk_p_hard_aj_spectra_bin[i]->Integral());
      h_trk_p_match_aj_spectra_bin[i]->Scale(auau_weights[i]/h_trk_p_match_aj_spectra_bin[i]->Integral());
      h_trk_m_hard_aj_spectra_bin[i]->Scale(auau_weights[i]/h_trk_m_hard_aj_spectra_bin[i]->Integral());
      h_trk_m_match_aj_spectra_bin[i]->Scale(auau_weights[i]/h_trk_m_match_aj_spectra_bin[i]->Integral());
    }
    
    std::vector<TH1D*> h_pp_hard_aj_spectra = AddBins(h_pp_hard_aj_spectra_bin, centrality_5);
    std::vector<TH1D*> h_pp_match_aj_spectra = AddBins(h_pp_match_aj_spectra_bin, centrality_5);
    
    std::vector<TH1D*> h_tow_p_hard_aj_spectra = AddBins(h_tow_p_hard_aj_spectra_bin, centrality_5);
    std::vector<TH1D*> h_tow_p_match_aj_spectra = AddBins(h_tow_p_match_aj_spectra_bin, centrality_5);
    std::vector<TH1D*> h_tow_m_hard_aj_spectra = AddBins(h_tow_m_hard_aj_spectra_bin, centrality_5);
    std::vector<TH1D*> h_tow_m_match_aj_spectra = AddBins(h_tow_m_match_aj_spectra_bin, centrality_5);
    
    std::vector<TH1D*> h_trk_p_hard_aj_spectra = AddBins(h_trk_p_hard_aj_spectra_bin, centrality_5);
    std::vector<TH1D*> h_trk_p_match_aj_spectra = AddBins(h_trk_p_match_aj_spectra_bin, centrality_5);
    std::vector<TH1D*> h_trk_m_hard_aj_spectra = AddBins(h_trk_m_hard_aj_spectra_bin, centrality_5);
    std::vector<TH1D*> h_trk_m_match_aj_spectra = AddBins(h_trk_m_match_aj_spectra_bin, centrality_5);
    
    std::vector<TH1D*> h_auau_hard_aj_spectra = SplitByCentrality(h_auau_hard_aj, centrality_5);
    std::vector<TH1D*> h_auau_match_aj_spectra = SplitByCentrality(h_auau_match_aj, centrality_5);
    
    // re-normalize
    for (int i = 0; i < h_auau_hard_aj_spectra.size(); ++i) {
      
      h_auau_hard_aj_spectra[i]->RebinX(2);
      h_auau_match_aj_spectra[i]->RebinX(2);
      
      h_auau_hard_aj_spectra[i]->Scale(1.0 / h_auau_hard_aj_spectra[i]->Integral());
      h_auau_match_aj_spectra[i]->Scale(1.0 / h_auau_match_aj_spectra[i]->Integral());
      
      h_pp_hard_aj_spectra[i]->Scale(1.0 / h_pp_hard_aj_spectra[i]->Integral());
      h_pp_match_aj_spectra[i]->Scale(1.0 / h_pp_match_aj_spectra[i]->Integral());
      
      h_tow_p_hard_aj_spectra[i]->Scale(1.0 / h_tow_p_hard_aj_spectra[i]->Integral());
      h_tow_p_match_aj_spectra[i]->Scale(1.0 / h_tow_p_match_aj_spectra[i]->Integral());
      h_tow_m_hard_aj_spectra[i]->Scale(1.0 / h_tow_m_hard_aj_spectra[i]->Integral());
      h_tow_m_match_aj_spectra[i]->Scale(1.0 / h_tow_m_match_aj_spectra[i]->Integral());
      
      h_trk_p_hard_aj_spectra[i]->Scale(1.0 / h_trk_p_hard_aj_spectra[i]->Integral());
      h_trk_p_match_aj_spectra[i]->Scale(1.0 / h_trk_p_match_aj_spectra[i]->Integral());
      h_trk_m_hard_aj_spectra[i]->Scale(1.0 / h_trk_m_hard_aj_spectra[i]->Integral());
      h_trk_m_match_aj_spectra[i]->Scale(1.0 / h_trk_m_match_aj_spectra[i]->Integral());
    }
    

    Overlay1D(h_auau_hard_aj_spectra, centrality_5_string, hOpts, cOpts, out_loc, "auau_hard_aj", "",
              "A_{J}", "fraction", "Centrality");
    Overlay1D(h_auau_match_aj_spectra, centrality_5_string, hOpts, cOpts, out_loc, "auau_match_aj", "",
              "A_{J}", "fraction", "Centrality");
    Overlay1D(h_pp_hard_aj_spectra, centrality_5_string, hOpts, cOpts, out_loc, "pp_hard_aj", "",
              "A_{J}", "fraction", "Centrality");
    Overlay1D(h_pp_match_aj_spectra, centrality_5_string, hOpts, cOpts, out_loc, "pp_match_aj", "",
              "A_{J}", "fraction", "Centrality");

    Overlay1D(h_tow_p_hard_aj_spectra, centrality_5_string, hOpts, cOpts, out_loc, "tow_p_hard_aj", "",
              "A_{J}", "fraction", "Centrality");
    Overlay1D(h_tow_p_match_aj_spectra, centrality_5_string, hOpts, cOpts, out_loc, "tow_p_match_aj", "",
              "A_{J}", "fraction", "Centrality");
    Overlay1D(h_tow_m_hard_aj_spectra, centrality_5_string, hOpts, cOpts, out_loc, "tow_m_hard_aj", "",
              "A_{J}", "fraction", "Centrality");
    Overlay1D(h_tow_m_match_aj_spectra, centrality_5_string, hOpts, cOpts, out_loc, "tow_m_match_aj", "",
              "A_{J}", "fraction", "Centrality");
    
    Overlay1D(h_trk_p_hard_aj_spectra, centrality_5_string, hOpts, cOpts, out_loc, "trk_p_hard_aj", "",
              "A_{J}", "fraction", "Centrality");
    Overlay1D(h_trk_p_match_aj_spectra, centrality_5_string, hOpts, cOpts, out_loc, "trk_p_match_aj", "",
              "A_{J}", "fraction", "Centrality");
    Overlay1D(h_trk_m_hard_aj_spectra, centrality_5_string, hOpts, cOpts, out_loc, "trk_m_hard_aj", "",
              "A_{J}", "fraction", "Centrality");
    Overlay1D(h_trk_m_match_aj_spectra, centrality_5_string, hOpts, cOpts, out_loc, "trk_m_match_aj", "",
              "A_{J}", "fraction", "Centrality");
    
    // add the containers to the dictionaries
    auau_hard_lead_pt_cent.insert({key, h_auau_hard_lead_pt_spectra});
    auau_hard_sub_pt_cent.insert({key, h_auau_hard_sub_pt_spectra});
    auau_match_lead_pt_cent.insert({key, h_auau_match_lead_pt_spectra});
    auau_match_sub_pt_cent.insert({key, h_auau_match_sub_pt_spectra});
    pp_hard_lead_pt_cent.insert({key, h_pp_hard_lead_pt_spectra});
    pp_hard_sub_pt_cent.insert({key, h_pp_hard_sub_pt_spectra});
    pp_match_lead_pt_cent.insert({key, h_pp_match_lead_pt_spectra});
    pp_match_sub_pt_cent.insert({key, h_pp_match_sub_pt_spectra});
    auau_hard_aj_cent.insert({key, h_auau_hard_aj_spectra});
    auau_match_aj_cent.insert({key, h_auau_match_aj_spectra});
    auau_hard_aj_cent.insert({key, h_pp_hard_aj_spectra});
    auau_match_aj_cent.insert({key, h_pp_match_aj_spectra});

    auau_hard_lead_rho_cent.insert({key, h_auau_hard_lead_rho_spectra});
    auau_hard_lead_sig_cent.insert({key, h_auau_hard_lead_sig_spectra});
    auau_hard_sub_rho_cent.insert({key, h_auau_hard_sub_rho_spectra});
    auau_hard_sub_sig_cent.insert({key, h_auau_hard_sub_sig_spectra});
    auau_match_lead_rho_cent.insert({key, h_auau_match_lead_rho_spectra});
    auau_match_lead_sig_cent.insert({key, h_auau_match_lead_sig_spectra});
    auau_match_sub_rho_cent.insert({key, h_auau_match_sub_rho_spectra});
    auau_match_sub_sig_cent.insert({key, h_auau_match_sub_sig_spectra});

    pp_hard_lead_rho_cent.insert({key, h_pp_hard_lead_rho_spectra});
    pp_hard_lead_sig_cent.insert({key, h_pp_hard_lead_sig_spectra});
    pp_hard_sub_rho_cent.insert({key, h_pp_hard_sub_rho_spectra});
    pp_hard_sub_sig_cent.insert({key, h_pp_hard_sub_sig_spectra});
    pp_match_lead_rho_cent.insert({key, h_pp_match_lead_rho_spectra});
    pp_match_lead_sig_cent.insert({key, h_pp_match_lead_sig_spectra});
    pp_match_sub_rho_cent.insert({key, h_pp_match_sub_rho_spectra});
    pp_match_sub_sig_cent.insert({key, h_pp_match_sub_sig_spectra});

    // now, printing some comparisons
    // make a folder for each centrality
    // make output directory
    for (int i = 0; i < h_auau_hard_aj_spectra.size(); ++i) {
      std::string out_loc_cent = out_loc + "/cent_" + std::to_string(i);
      boost::filesystem::path dir(out_loc_cent.c_str());
      boost::filesystem::create_directories(dir);
      
      // first overlay all the systematics
      Overlay1D(h_tow_p_hard_aj_spectra[i], h_tow_m_hard_aj_spectra[i], "+2% tower", "-2% tower",
                hOpts, cOpts, out_loc_cent, "aj_hard_tow_sys", "", "A_{J}", "fraction", "A_{J}");
      Overlay1D(h_tow_p_match_aj_spectra[i], h_tow_m_match_aj_spectra[i], "+2% tower", "-2% tower",
                hOpts, cOpts, out_loc_cent, "aj_match_tow_sys", "", "A_{J}", "fraction", "A_{J}");
      Overlay1D(h_trk_p_hard_aj_spectra[i], h_trk_m_hard_aj_spectra[i], "+ track eff", "- track eff",
                hOpts, cOpts, out_loc_cent, "aj_hard_trk_sys", "", "A_{J}", "fraction", "A_{J}");
      Overlay1D(h_trk_p_match_aj_spectra[i], h_trk_m_match_aj_spectra[i], "+ track eff", "- track eff",
                hOpts, cOpts, out_loc_cent, "aj_match_trk_sys", "", "A_{J}", "fraction", "A_{J}");
      
      
      TGraphErrors* err_hard = GetSystematic(h_pp_hard_aj_spectra[i], h_tow_p_hard_aj_spectra[i], h_tow_m_hard_aj_spectra[i],
                                                 h_trk_p_hard_aj_spectra[i], h_trk_m_hard_aj_spectra[i]);
      TGraphErrors* err_match = GetSystematic(h_pp_match_aj_spectra[i], h_tow_p_match_aj_spectra[i], h_tow_m_match_aj_spectra[i],
                                                  h_trk_p_match_aj_spectra[i], h_trk_m_match_aj_spectra[i]);

      Overlay1D(h_auau_hard_aj_spectra[i], h_pp_hard_aj_spectra[i], 0.0, 0.25, 0.0, 0.9, "AuAu hard A_{J}", "PP hard A_{J}",
                 hOpts, cOpts, out_loc_cent, "aj_hard", "", "A_{J}", "fraction", "A_{J}");
      Overlay1D(h_auau_match_aj_spectra[i], h_pp_match_aj_spectra[i], 0.0, 0.23, 0.0, 0.9, "AuAu matched A_{J}", "PP matched A_{J}",
                 hOpts, cOpts, out_loc_cent, "aj_match", "", "A_{J}", "fraction", "A_{J}");
      
      Overlay1D(h_auau_hard_aj_spectra[i], h_pp_hard_aj_spectra[i], err_hard, 0.0, 0.25, 0.0, 0.9, "AuAu hard A_{J}", "PP hard A_{J}",
                "systematics", hOpts, cOpts, out_loc_cent, "aj_hard", "", "A_{J}", "fraction", "A_{J}");
      Overlay1D(h_auau_match_aj_spectra[i], h_pp_match_aj_spectra[i], err_match, 0.0, 0.3, 0.0, 0.9, "AuAu matched A_{J}", "PP matched A_{J}",
                "systematics", hOpts, cOpts, out_loc_cent, "aj_match", "", "A_{J}", "fraction", "A_{J}");
    }
    
    string out_loc_aj = out_loc + "/aj_plots";
    boost::filesystem::path aj_dir(out_loc_aj);
    boost::filesystem::create_directories(aj_dir);
    
    TFile* kolja = new TFile("kolja.root", "READ");
    TH1D* kolja_au = (TH1D*) kolja->Get("AuAuAJ_hi");
    hOpts.SetHistogram(kolja_au);
    kolja_au->SetMarkerStyle(29);
    kolja_au->SetMarkerSize(2);
    kolja_au->SetMarkerColor(kRed);
    kolja_au->SetLineColor(kRed);
    kolja_au->GetXaxis()->SetTitle("A_{J}");
    kolja_au->GetYaxis()->SetTitle("fraction");
    
    TH1D* kolja_pp = (TH1D*) kolja->Get("ppInAuAuAJ_hi");
    kolja_pp->SetMarkerStyle(30);
    kolja_pp->SetMarkerSize(2);
    
    
    string name_aj1 = out_loc + "/aj_plots/aj_0_20.pdf";
    string name_aj1G = out_loc + "/aj_plots/aj_0_20.gif";
    TCanvas* c1 = new TCanvas();
    cOpts.SetMargins(c1);
    
    kolja_au->GetXaxis()->SetRangeUser(0.0, 0.7);
    kolja_au->Draw();
    //kolja_pp->Draw("same");
    
    h_auau_hard_aj_spectra[0]->SetMarkerSize(1.7);
    h_auau_hard_aj_spectra[0]->SetLineColor(kBlue);
    h_auau_hard_aj_spectra[0]->SetMarkerColor(kBlue);
    h_auau_hard_aj_spectra[0]->GetYaxis()->SetRangeUser(0.0, 0.25);
    h_auau_hard_aj_spectra[0]->GetXaxis()->SetTitle("A_{J}");
    h_auau_hard_aj_spectra[0]->GetYaxis()->SetTitle("fraction");
    h_auau_hard_aj_spectra[0]->Draw("SAME");
    
    TLegend* leg1 = new TLegend(0.55, 0.7, 0.9, 0.88);
    leg1->SetTextSize(0.032);
    leg1->AddEntry(h_auau_hard_aj_spectra[0], "Run 14 Au+Au", "lep");
    leg1->AddEntry(kolja_au, "Run 7 Au+Au", "lep");
    //leg1->AddEntry(kolja_pp, "Run 6 p+p HT + Au+Au 0-20%", "lep");
    leg1->Draw();
    
    TPaveText *t = new TPaveText(0.65, 0.5, 0.88, 0.65, "NB NDC");
    t->SetFillStyle(0);
    t->SetBorderSize(0);
    string cent_string_1 = centrality_5_string[0] + " Centrality";
    t->AddText("p_{T}^{const} > 2.0 GeV/c");
    t->AddText("p+p eff. corrected to Run 7 Au+Au 0-20%" );
    //t->Draw();
    
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.045);
    latex.SetTextColor(kGray+3);
    latex.SetTextColor(kRed+3);
    latex.DrawLatex( 0.15, 0.85, "STAR Preliminary");
    
    TLatex latex1;
    latex1.SetNDC();
    latex1.SetTextSize(0.04);
    latex1.SetTextColor(kBlack);
    latex1.DrawLatex( 0.67, 0.44, cent_string_1.c_str());
    
    TLatex latex2;
    latex2.SetNDC();
    latex2.SetTextSize(0.04);
    latex2.SetTextColor(kBlack);
    latex2.DrawLatex( 0.25, 0.42, "p_{T}^{const} > 2.0 GeV/c");
    latex2.DrawLatex( 0.25, 0.36, "p_{T}^{lead} > 20.0 GeV/c");
    latex2.DrawLatex( 0.25, 0.3, "p_{T}^{sublead} > 10.0 GeV/c");
    //latex2.DrawLatex( 0.22, 0.24, "p+p eff. corrected to Run 7");
    //latex2.DrawLatex( 0.29, 0.2, "Au+Au 0-20%" );
    
    c1->SaveAs(name_aj1.c_str());
    c1->SaveAs(name_aj1G.c_str());
    
//    TCanvas *c2 = new TCanvas();
//    cOpts.SetMargins(c2);
//
//    h_auau_hard_aj_spectra[1]->GetXaxis()->SetRangeUser(0.0, 0.7);
//    h_auau_hard_aj_spectra[1]->SetMarkerStyle(20);
//    h_auau_hard_aj_spectra[1]->SetMarkerSize(1.7);
//    h_auau_hard_aj_spectra[1]->SetLineColor(kMagenta);
//    h_auau_hard_aj_spectra[1]->SetMarkerColor(kMagenta);
//    h_auau_hard_aj_spectra[1]->GetYaxis()->SetRangeUser(0.0, 0.25);
//    h_auau_hard_aj_spectra[1]->GetXaxis()->SetTitle("A_{J}");
//    h_auau_hard_aj_spectra[1]->GetYaxis()->SetTitle("fraction");
//    h_auau_hard_aj_spectra[1]->Draw();
//
//    TLegend* leg2 = new TLegend(0.65, 0.8, 0.88, 0.88);
//    leg2->AddEntry(h_auau_hard_aj_spectra[1], "Run 14 Au+Au", "lep");
//    leg2->Draw();
//
//    latex.DrawLatex( 0.15, 0.85, "STAR Preliminary");
//    string cent_string_2 = centrality_5_string[1] + " Centrality";
//    latex1.DrawLatex( 0.67, 0.44, cent_string_2.c_str());
//
//    string name_aj2 = out_loc + "/aj_plots/aj_20_40.pdf";
//    c2->SaveAs(name_aj2.c_str());
//
//
//    TCanvas *c3 = new TCanvas();
//    cOpts.SetMargins(c3);
//
//    h_auau_hard_aj_spectra[2]->GetXaxis()->SetRangeUser(0.0, 0.7);
//    h_auau_hard_aj_spectra[2]->SetMarkerStyle(20);
//    h_auau_hard_aj_spectra[2]->SetMarkerSize(1.7);
//    h_auau_hard_aj_spectra[2]->SetLineColor(kGreen+1);
//    h_auau_hard_aj_spectra[2]->SetMarkerColor(kGreen+1);
//    h_auau_hard_aj_spectra[2]->GetYaxis()->SetRangeUser(0.0, 0.25);
//    h_auau_hard_aj_spectra[2]->GetXaxis()->SetTitle("A_{J}");
//    h_auau_hard_aj_spectra[2]->GetYaxis()->SetTitle("fraction");
//    h_auau_hard_aj_spectra[2]->Draw();
//
//    TLegend* leg3 = new TLegend(0.65, 0.8, 0.88, 0.88);
//    leg3->AddEntry(h_auau_hard_aj_spectra[2], "Run 14 Au+Au", "lep");
//    leg3->Draw();
//
//    latex.DrawLatex( 0.15, 0.85, "STAR Preliminary");
//    string cent_string_3 = centrality_5_string[2] + " Centrality";
//    latex1.DrawLatex( 0.67, 0.44, cent_string_3.c_str());
//
//    string name_aj3 = out_loc + "/aj_plots/aj_50_70.pdf";
//    c3->SaveAs(name_aj3.c_str());
    
    
    TCanvas *c2 = new TCanvas();
    cOpts.SetMargins(c2);
    
    h_auau_hard_aj_spectra[0]->GetXaxis()->SetRangeUser(0.0, 0.7);
    h_auau_hard_aj_spectra[0]->SetMarkerSize(1.7);
    h_auau_hard_aj_spectra[0]->SetLineColor(kBlue);
    h_auau_hard_aj_spectra[0]->SetMarkerColor(kBlue);
    h_auau_hard_aj_spectra[0]->GetYaxis()->SetRangeUser(0.0, 0.25);
    h_auau_hard_aj_spectra[0]->GetXaxis()->SetTitle("A_{J}");
    h_auau_hard_aj_spectra[0]->GetYaxis()->SetTitle("fraction");
    h_auau_hard_aj_spectra[0]->Draw("");
    std::cout << "entries in 0-20: " << h_auau_hard_aj_spectra[0]->GetEntries() << std::endl;
    
    h_auau_hard_aj_spectra[1]->GetXaxis()->SetRangeUser(0.0, 0.7);
    h_auau_hard_aj_spectra[1]->SetMarkerStyle(20);
    h_auau_hard_aj_spectra[1]->SetMarkerSize(1.7);
    h_auau_hard_aj_spectra[1]->SetLineColor(kMagenta);
    h_auau_hard_aj_spectra[1]->SetMarkerColor(kMagenta);
    h_auau_hard_aj_spectra[1]->GetYaxis()->SetRangeUser(0.0, 0.25);
    h_auau_hard_aj_spectra[1]->GetXaxis()->SetTitle("A_{J}");
    h_auau_hard_aj_spectra[1]->GetYaxis()->SetTitle("fraction");
    h_auau_hard_aj_spectra[1]->Draw("SAME");
    std::cout << "entries in 20-40: " << h_auau_hard_aj_spectra[1]->GetEntries() << std::endl;
    
    h_auau_hard_aj_spectra[2]->GetXaxis()->SetRangeUser(0.0, 0.7);
    h_auau_hard_aj_spectra[2]->SetMarkerStyle(20);
    h_auau_hard_aj_spectra[2]->SetMarkerSize(1.7);
    h_auau_hard_aj_spectra[2]->SetLineColor(kGreen+1);
    h_auau_hard_aj_spectra[2]->SetMarkerColor(kGreen+1);
    h_auau_hard_aj_spectra[2]->GetYaxis()->SetRangeUser(0.0, 0.25);
    h_auau_hard_aj_spectra[2]->GetXaxis()->SetTitle("A_{J}");
    h_auau_hard_aj_spectra[2]->GetYaxis()->SetTitle("fraction");
    h_auau_hard_aj_spectra[2]->Draw("SAME");
    std::cout << "entries in 50-70: " << h_auau_hard_aj_spectra[2]->GetEntries() << std::endl;
    
    TLegend* leg2 = new TLegend(0.55, 0.65, 0.88, 0.88);
    leg2->SetTextSize(0.032);
    leg2->AddEntry(h_auau_hard_aj_spectra[0], "Run 14 0-20%", "lep");
    leg2->AddEntry(h_auau_hard_aj_spectra[1], "Run 14 20-40%", "lep");
    leg2->AddEntry(h_auau_hard_aj_spectra[2], "Run 14 50-70%", "lep");
    leg2->Draw();
    
    latex.DrawLatex( 0.15, 0.85, "STAR Preliminary");
    
    latex2.SetNDC();
    latex2.SetTextSize(0.04);
    latex2.SetTextColor(kBlack);
    latex2.DrawLatex( 0.25, 0.42, "p_{T}^{const} > 2.0 GeV/c");
    latex2.DrawLatex( 0.25, 0.36, "p_{T}^{lead} > 20.0 GeV/c");
    latex2.DrawLatex( 0.25, 0.3, "p_{T}^{sublead} > 10.0 GeV/c");
    
    string name_aj2 = out_loc + "/aj_plots/aj_run14.pdf";
    string name_aj2G = out_loc + "/aj_plots/aj_run14.gif";
    c2->SaveAs(name_aj2.c_str());
    c2->SaveAs(name_aj2G.c_str());
    
  }
  
  return 0;
}

