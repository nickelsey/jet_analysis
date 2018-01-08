// print_full_results.cc

#include "jet_analysis/util/root_print_routines.hh"
#include "jet_analysis/util/arg_helper.hh"

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

std::vector<TH1D*> SplitByRefMult(TH2D* h, std::vector<int>& refmult) {
  
  std::vector<TH1D*> ret;
  int bound = h->GetXaxis()->GetNbins();
  for (int i = 0; i < refmult.size(); ++i) {
    string name = string(h->GetName()) + std::to_string(i);
    TH1D* tmp = h->ProjectionY(name.c_str(),
                               h->GetXaxis()->FindBin(refmult[i]),
                               bound);
    bound = h->GetXaxis()->FindBin(refmult[i]) - 1;
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

// command line option results
struct Options {
  string auau_file  = "";
  string pp_file    = "";
  string out_loc    = "results";
  string out_prefix = "";
  int cent_low      = 0;
  int cent_high     = 8;
};

int main(int argc, char* argv[]) {
  
  // set drawing preferences for histograms and graphs
  gStyle->SetOptStat(false);
  gStyle->SetOptFit(false);
  gStyle->SetOptTitle(1);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetHatchesSpacing(1.0);
  gStyle->SetHatchesLineWidth(2);
  
  Options opts;
  for (int i = 1; i < argc; ++i) {
    if (ParseStrFlag(string(argv[i]), "--auau", &opts.auau_file) ||
        ParseStrFlag(string(argv[i]), "--pp", &opts.pp_file) ||
        ParseStrFlag(string(argv[i]), "--outputDir", &opts.out_loc) ||
        ParseStrFlag(string(argv[i]), "--filePrefix", &opts.out_prefix) ||
        ParseIntFlag(string(argv[i]), "--centLow", &opts.cent_low) ||
        ParseIntFlag(string(argv[i]), "--centHigh", &opts.cent_high)) continue;
    std::cerr << "Unknown command line option: " << argv[i] << std::endl;
    return 1;
  }
  
  // check to make sure we have valid inputs
  if (!boost::filesystem::exists(opts.pp_file)) {
    std::cout << "AuAu input file " << opts.auau_file;
    std::cout << " doesn't exist: exiting" << std::endl;
    return 1;
  }
  if (!boost::filesystem::exists(opts.pp_file)) {
    std::cout << "PP input file " << opts.pp_file;
    std::cout << " doesn't exist: exiting" << std::endl;
    return 1;
  }
  
  // read in the two files
  TFile auau_file(opts.auau_file.c_str(), "READ");
  TFile pp_file(opts.pp_file.c_str(), "READ");
  
  // now we'll get the trees from the files, ignoring any objects
  // in the file that don't conform to the naming conventions from
  // the DijetWorker
  std::vector<string> keys;
  std::unordered_map<std::string, TTree*> auau_trees;
  std::unordered_map<std::string, TTree*> pp_trees;
  
  GetTreesFromFile(auau_file, auau_trees);
  GetTreesFromFile(pp_file, pp_trees);
  
  // match keys
  for (auto entry : auau_trees) {
    if (pp_trees.find(entry.first) == pp_trees.end())
      auau_trees.erase(entry.first);
  }
  for (auto entry : pp_trees) {
    if (auau_trees.find(entry.first) == auau_trees.end())
      pp_trees.erase(entry.first);
    else
      keys.push_back(entry.first);
  }
  
  // for now, use hard coded centrality definition for run 14
  std::vector<int> refcent_def{420, 364, 276, 212, 156, 108, 68, 44, 28, 12, 0};
  std::vector<string> refcent_string{"0-5%", "5-10%", "10-20%", "20-30%",
                                     "30-40%", "40-50%", "40-60%", "60-70%",
                                     "70-80%", "80-90%", "90-100%"};
  std::vector<string> refcent_alt_string{"0-5%", "5-10%", "10-20%", "20-30%",
                                         "30-40%", "40-50%", "40-60%", "60-70%",
                                         "70-80%"};
  
  // loop over all matched trees
  for (auto key : keys) {
    // make output directory
    string out_loc = opts.out_loc += "/" + key;
    boost::filesystem::path dir(out_loc.c_str());
    boost::filesystem::create_directories(dir);
    
    // and create the file name prefix
    string file_prefix = out_loc + "/" + opts.out_prefix;
    
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
    TTreeReaderValue<int> embed_eventid;
    TTreeReaderValue<int> embed_runid;
    TTreeReaderValue<int> embed_refmult;
    TTreeReaderValue<int> embed_grefmult;
    TTreeReaderValue<double> embed_refmultcorr;
    TTreeReaderValue<double> embed_grefmultcorr;
    TTreeReaderValue<int> embed_cent;
    TTreeReaderValue<double> embed_rp;
    TTreeReaderValue<double> embed_zdcrate;
    TTreeReaderValue<double> embed_vz;
    bool pp_embedded = false;
    if (pp_tree->GetBranch("embed_eventid")) {
      pp_embedded = true;
      embed_eventid = TTreeReaderValue<int>(pp_reader, "embed_eventid");
      embed_runid = TTreeReaderValue<int>(pp_reader, "embed_runid");
      embed_refmult = TTreeReaderValue<int>(pp_reader, "embed_refmult");
      embed_grefmult = TTreeReaderValue<int>(pp_reader, "embed_grefmult");
      embed_refmultcorr = TTreeReaderValue<double>(pp_reader, "embed_refmultcorr");
      embed_grefmultcorr = TTreeReaderValue<double>(pp_reader, "embed_grefmultcorr");
      embed_cent = TTreeReaderValue<int>(pp_reader, "embed_cent");
      embed_rp = TTreeReaderValue<double>(pp_reader, "embed_rp");
      embed_zdcrate = TTreeReaderValue<double>(pp_reader, "embed_zdcrate");
      embed_vz = TTreeReaderValue<double>(pp_reader, "embed_vz");
    }
    
    // -----------------
    // define histograms
    // -----------------
    
    // jet pt
    TH2D* auau_hard_lead_pt = new TH2D("auauhardleadpt", "p_{T}", 800, 0.5, 800.5, 100, 0, 100);
    TH2D* auau_hard_sub_pt = new TH2D("auauhardsubpt", "p_{T}", 800, 0.5, 800.5, 100, 0, 100);
    TH2D* auau_match_lead_pt = new TH2D("auaumatchleadpt", "p_{T}", 800, 0.5, 800.5, 100, 0, 100);
    TH2D* auau_match_sub_pt = new TH2D("auaumatchsubpt", "p_{T}", 800, 0.5, 800.5, 100, 0, 100);
    TH2D* pp_hard_lead_pt = new TH2D("pphardleadpt", "p_{T}", 800, 0.5, 800.5, 100, 0, 100);
    TH2D* pp_hard_sub_pt = new TH2D("pphardsubpt", "p_{T}", 800, 0.5, 800.5, 100, 0, 100);
    TH2D* pp_match_lead_pt = new TH2D("ppmatchleadpt", "p_{T}", 800, 0.5, 800.5, 100, 0, 100);
    TH2D* pp_match_sub_pt = new TH2D("ppmatchsubpt", "p_{T}", 800, 0.5, 800.5, 100, 0, 100);
    
    // aj
    TH2D* auau_hard_aj = new TH2D("auauhardaj", "A_{J}", 800, 0.5, 800.5, 30, 0, 0.9);
    TH2D* auau_match_aj = new TH2D("auaumatchaj", "A_{J}", 800, 0.5, 800.5, 30, 0, 0.9);
    TH2D* pp_hard_aj = new TH2D("pphardaj", "A_{J}", 800, 0.5, 800.5, 30, 0, 0.9);
    TH2D* pp_match_aj = new TH2D("ppmatchaj", "A_{J}", 800, 0.5, 800.5, 30, 0, 0.9);
    
    // dphi
    TH2D* auau_dphi = new TH2D("auaudphi", "d#phi", 800, 0.5, 800.5, 100, 0, 2*TMath::Pi());
    TH2D* pp_dphi = new TH2D("ppdphi", "d#phi", 800, 0.5, 800.5, 100, 0, 2*TMath::Pi());
    
    // lead jet - reaction plane
    TH2D* auau_hard_lead_rp = new TH2D("auauhardleaddphi", "phi - rp", 9, -0.5, 8.5, 100, 0, 2*TMath::Pi());
    
    
    // loop over the data & fill histograms
    while (auau_reader.Next()) {
      
      // auau jet pt
      auau_hard_lead_pt->Fill(*auau_refmult, (*auau_jl).Pt());
      auau_hard_sub_pt->Fill(*auau_refmult, (*auau_js).Pt());
      auau_match_lead_pt->Fill(*auau_refmult,(*auau_jlm).Pt());
      auau_match_sub_pt->Fill(*auau_refmult, (*auau_jsm).Pt());
      
      // auau Aj
      auau_hard_aj->Fill(*auau_refmult,
                         fabs((*auau_jl).Pt() - (*auau_js).Pt())/((*auau_jl).Pt() + (*auau_js).Pt()));
      auau_match_aj->Fill(*auau_refmult,
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
      
      auau_dphi->Fill(*auau_refmult, dphi);
      auau_hard_lead_rp->Fill(*auau_cent, dphi_rp);
      
    }
    
    while (pp_reader.Next()) {
      if (pp_embedded) {
        
        // pp single jet pt
        pp_hard_lead_pt->Fill(*pp_refmult + *embed_refmult, (*pp_jl).Pt());
        pp_hard_sub_pt->Fill(*pp_refmult + *embed_refmult, (*pp_js).Pt());
        pp_match_lead_pt->Fill(*pp_refmult + *embed_refmult, (*pp_jlm).Pt());
        pp_match_sub_pt->Fill(*pp_refmult + *embed_refmult, (*pp_jsm).Pt());
        
        // pp Aj
        pp_hard_aj->Fill(*pp_refmult + *embed_refmult,
                         fabs((*pp_jl).Pt() - (*pp_js).Pt())/((*pp_jl).Pt() + (*pp_js).Pt()));
        pp_match_aj->Fill(*pp_refmult + *embed_refmult,
                          fabs((*pp_jlm).Pt() - (*pp_jsm).Pt())/((*pp_jlm).Pt() + (*pp_jsm).Pt()));
        
        // pp dphi
        double dphi = (*pp_jl).Phi() - (*pp_js).Phi();
        //rotate dphi to be in [0,2*pi]
        while (dphi < 0)
          dphi += 2 * TMath::Pi();
        while (dphi > 2.0 * TMath::Pi())
          dphi -= 2.0 * TMath::Pi();
        
        pp_dphi->Fill(*pp_refmult + *embed_refmult, dphi);
      }
      else {
        
        // pp single jet pt
        pp_hard_lead_pt->Fill(*pp_refmult, (*pp_jl).Pt());
        pp_hard_sub_pt->Fill(*pp_refmult, (*pp_js).Pt());
        pp_match_lead_pt->Fill(*pp_refmult, (*pp_jlm).Pt());
        pp_match_sub_pt->Fill(*pp_refmult, (*pp_jsm).Pt());
        
        // pp Aj
        pp_hard_aj->Fill(*pp_refmult,
                         ((*pp_jl).Pt() - (*pp_js).Pt())/((*pp_jl).Pt() + (*pp_js).Pt()));
        pp_match_aj->Fill(*pp_refmult,
                          ((*pp_jlm).Pt() - (*pp_jsm).Pt())/((*pp_jlm).Pt() + (*pp_jsm).Pt()));
        
        // pp dphi
        double dphi = (*pp_jl).Phi() - (*pp_js).Phi();
        //rotate dphi to be in [0,2*pi]
        while (dphi < 0)
        dphi += 2 * TMath::Pi();
        while (dphi > 2.0 * TMath::Pi())
        dphi -= 2.0 * TMath::Pi();
        
        pp_dphi->Fill(*pp_refmult, dphi);
      }
    }
    
    // print dphi
    auau_dphi->Scale(1.0/auau_dphi->Integral());
    pp_dphi->Scale(1.0/pp_dphi->Integral());
    
    Overlay1D((TH1D*)auau_dphi->ProjectionY(), (TH1D*)pp_dphi->ProjectionY(),
              "Au+Au d#phi lead-sub", "P+P d#phi lead-sub",
              out_loc, "auau_pp_dphi", "", "d#phi", "fraction", false, false, false, "");
    Print2DSimple(auau_hard_lead_rp, out_loc, "auau_dphi_rp", "AuAu d#phi lead jet - rp",
                  "Centrality", "d#phi(lead-rp)", false, false, false);
    std::vector<TH1D*> rp_by_cent = SplitByBin(auau_hard_lead_rp);
    for (int i = 0; i < rp_by_cent.size(); ++i) {
      rp_by_cent[i]->Scale(1.0/rp_by_cent[i]->Integral());
      rp_by_cent[i]->RebinX(4);
    }
    Overlay1D(rp_by_cent, refcent_alt_string, out_loc, "auau_cent_rp", "",
              "d#phi", "fraction", false, false, true, "Centrality");
    Overlay1D(rp_by_cent[0], rp_by_cent[1], "0-5%", "5-10%", out_loc, "auau_cent_rp_restricted", "",
              "d#phi", "fraction", false, false, true, "Centrality");
    
    // extract pt spectra in centrality bins
    std::vector<TH1D*> auau_hard_lead_pt_spectra = SplitByRefMult(auau_hard_lead_pt,
                                                                  refcent_def);
    std::vector<TH1D*> auau_hard_sub_pt_spectra = SplitByRefMult(auau_hard_sub_pt,
                                                                 refcent_def);
    std::vector<TH1D*> auau_match_lead_pt_spectra = SplitByRefMult(auau_match_lead_pt,
                                                                   refcent_def);
    std::vector<TH1D*> auau_match_sub_pt_spectra = SplitByRefMult(auau_match_sub_pt,
                                                                  refcent_def);
    std::vector<TH1D*> pp_hard_lead_pt_spectra = SplitByRefMult(pp_hard_lead_pt,
                                                                refcent_def);
    std::vector<TH1D*> pp_hard_sub_pt_spectra = SplitByRefMult(pp_hard_sub_pt,
                                                               refcent_def);
    std::vector<TH1D*> pp_match_lead_pt_spectra = SplitByRefMult(pp_match_lead_pt,
                                                                 refcent_def);
    std::vector<TH1D*> pp_match_sub_pt_spectra = SplitByRefMult(pp_match_sub_pt,
                                                                refcent_def);
    
    // normalize pt spectra
    for (int i = 0; i < auau_hard_lead_pt_spectra.size(); ++i) {
      auau_hard_lead_pt_spectra[i]->Scale(1.0/auau_hard_lead_pt_spectra[i]->Integral());
      auau_hard_sub_pt_spectra[i]->Scale(1.0/auau_hard_sub_pt_spectra[i]->Integral());
      auau_match_lead_pt_spectra[i]->Scale(1.0/auau_match_lead_pt_spectra[i]->Integral());
      auau_match_sub_pt_spectra[i]->Scale(1.0/auau_match_sub_pt_spectra[i]->Integral());
      pp_hard_lead_pt_spectra[i]->Scale(1.0/pp_hard_lead_pt_spectra[i]->Integral());
      pp_hard_sub_pt_spectra[i]->Scale(1.0/pp_hard_sub_pt_spectra[i]->Integral());
      pp_match_lead_pt_spectra[i]->Scale(1.0/pp_match_lead_pt_spectra[i]->Integral());
      pp_match_sub_pt_spectra[i]->Scale(1.0/pp_match_sub_pt_spectra[i]->Integral());
    }
    
    // print pt spectra
    Overlay1D(auau_hard_lead_pt_spectra, refcent_string, out_loc,
              "auau_hard_lead_pt", "", "p_{T}", "fraction", false, true, true, "Centrality");
    Overlay1D(auau_hard_sub_pt_spectra, refcent_string, out_loc,
              "auau_hard_sub_pt", "", "p_{T}", "fraction", false, true, true, "Centrality");
    Overlay1D(auau_match_lead_pt_spectra, refcent_string, out_loc,
              "auau_match_lead_pt", "", "p_{T}", "fraction", false, true, true, "Centrality");
    Overlay1D(auau_match_sub_pt_spectra, refcent_string, out_loc,
              "auau_match_sub_pt", "", "p_{T}", "fraction", false, true, true, "Centrality");
    Overlay1D(pp_hard_lead_pt_spectra, refcent_string, out_loc,
              "pp_hard_lead_pt", "", "p_{T}", "fraction", false, true, true, "Centrality");
    Overlay1D(pp_hard_sub_pt_spectra, refcent_string, out_loc,
              "pp_hard_sub_pt", "", "p_{T}", "fraction", false, true, true, "Centrality");
    Overlay1D(pp_match_lead_pt_spectra, refcent_string, out_loc,
              "pp_match_lead_pt", "", "p_{T}", "fraction", false, true, true, "Centrality");
    Overlay1D(pp_match_sub_pt_spectra, refcent_string, out_loc,
              "pp_match_sub_pt", "", "p_{T}", "fraction", false, true, true, "Centrality");
    
    // for Aj
    std::vector<TH1D*> auau_hard_aj_spectra = SplitByRefMult(auau_hard_aj, refcent_def);
    std::vector<TH1D*> auau_match_aj_spectra = SplitByRefMult(auau_match_aj, refcent_def);
    std::vector<TH1D*> pp_hard_aj_spectra = SplitByRefMult(pp_hard_aj, refcent_def);
    std::vector<TH1D*> pp_match_aj_spectra = SplitByRefMult(pp_match_aj, refcent_def);
    
    // normalize
    for (int i = 0; i < auau_hard_aj_spectra.size(); ++i) {
      auau_hard_aj_spectra[i]->RebinX(2);
      auau_match_aj_spectra[i]->RebinX(2);
      pp_hard_aj_spectra[i]->RebinX(2);
      pp_match_aj_spectra[i]->RebinX(2);
      auau_hard_aj_spectra[i]->Scale(1.0/auau_hard_aj_spectra[i]->Integral());
      auau_match_aj_spectra[i]->Scale(1.0/auau_match_aj_spectra[i]->Integral());
      pp_hard_aj_spectra[i]->Scale(1.0/pp_hard_aj_spectra[i]->Integral());
      pp_match_aj_spectra[i]->Scale(1.0/pp_match_aj_spectra[i]->Integral());
    }
    Overlay1D(auau_hard_aj_spectra, refcent_string, out_loc, "auau_hard_aj", "",
              "A_{J}", "fraction", false, false, true, "Centrality", 0.22, 0.0);
    Overlay1D(auau_match_aj_spectra, refcent_string, out_loc, "auau_match_aj", "",
              "A_{J}", "fraction", false, false, true, "Centrality", 0.22, 0.0);
    Overlay1D(pp_hard_aj_spectra, refcent_string, out_loc, "pp_hard_aj", "",
              "A_{J}", "fraction", false, false, true, "Centrality", 0.22, 0.0);
    Overlay1D(pp_match_aj_spectra, refcent_string, out_loc, "pp_match_aj", "",
              "A_{J}", "fraction", false, false, true, "Centrality", 0.22, 0.0);
    
  }
  
  return 0;
}
