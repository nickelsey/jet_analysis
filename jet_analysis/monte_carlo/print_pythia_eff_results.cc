// print_pythia_eff_results.cc

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "jet_analysis/util/string_util.hh"
#include "jet_analysis/util/root_print_routines.hh"
#include "jet_analysis/util/arg_helper.hh"

#include <iostream>
#include <vector>
#include <string>

using std::string;
struct Options {
  string in_file  = "tmp/pythia.root";    /* input root file */
  string out_dir  = "tmp";   /* output location */
  int nbins       = 6;       /* number of efficiency bins */
};

int main(int argc, char* argv[]) {
  
  std::vector<string> eff_strings{"100%", "90%", "80%", "70%", "60%",
                                  "50%", "40%", "30%", "20%", "10%"};
  
  // create histogram options
  histogramOpts hOpts;
  canvasOpts cOpts;
  
  // Histograms will calculate gaussian errors
  // -----------------------------------------
  TH1::SetDefaultSumw2( );
  TH2::SetDefaultSumw2( );
  TH3::SetDefaultSumw2( );
  
  // set drawing preferences for histograms and graphs
  gStyle->SetOptStat(false);
  gStyle->SetOptFit(false);
  gStyle->SetOptTitle(1);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetHatchesSpacing(1.0);
  gStyle->SetHatchesLineWidth(2);
  
  Options opts;
  for (int i = 1; i < argc; ++i) {
    if (ParseIntFlag(string(argv[i]), "--effBins", &opts.nbins) ||
        ParseStrFlag(string(argv[i]), "--outDir", &opts.out_dir) ||
        ParseStrFlag(string(argv[i]), "--inFile", &opts.in_file)) continue;
    std::cerr << "unknown command line option: " << argv[i] << std::endl;
    return 1;
  }
  
  // build output directory if it doesn't exist, using boost::filesystem
  if (opts.out_dir.empty())
    opts.out_dir = "tmp";
  boost::filesystem::path dir(opts.out_dir.c_str());
  boost::filesystem::create_directories(dir);
  
  TFile in(opts.in_file.c_str(), "READ");
  
  std::vector<TH1D*> h_aj;
  std::vector<TH1D*> h_count1;
  std::vector<TH1D*> h_count2;
  
  for (int i = 0; i < opts.nbins; ++i) {
    h_aj.push_back((TH1D*) in.Get(MakeString("aj_", i).c_str()));
    h_aj[i]->Scale(1.0/h_aj[i]->Integral());
    h_aj[i]->Rebin(2);
    h_count1.push_back((TH1D*) in.Get(MakeString("npa_", i).c_str()));
    h_count1[i]->Scale(1.0/h_count1[i]->Integral());
    h_count1[i]->Rebin(2);
    h_count2.push_back((TH1D*) in.Get(MakeString("npb_", i).c_str()));
    h_count2[i]->Scale(1.0/h_count2[i]->Integral());
    h_count2[i]->Rebin(2);
  }
  
  eff_strings.resize(opts.nbins);
  
  
  
  Overlay1D(h_aj, eff_strings,  hOpts, cOpts, opts.out_dir,
            "aj", "", "A_{J}", "fraction", "efficiency");
  Overlay1D(h_count1, eff_strings, hOpts, cOpts, opts.out_dir,
            "ncount_after", "", "nPart", "fraction", "nPart after efficiency");
  Overlay1D(h_count2, eff_strings, hOpts, cOpts, opts.out_dir,
            "ncount_before", "", "nPart", "fraction", "nPart before efficiency");
  return 0;
}
