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

using std::string;

struct Options {
  string input_y7    = "";       /* root file for y7 data*/
  string input_y14   = "";       /* root file for y14 data */
  string out_dir     = "tmp";    /* directory to save output in */
};

int main(int argc, char* argv[]) {
  
  // parse command line options
  // --------------------------
  Options opts;
  for (int i = 1; i < argc; ++i) {
    if (ParseStrFlag(string(argv[i]), "--inputY7", &opts.input_y7) ||
        ParseStrFlag(string(argv[i]), "--inputY14", &opts.input_y14) ||
        ParseStrFlag(string(argv[i]), "--outDir", &opts.out_dir)) continue;
    std::cerr << "Unknown command line option: " << argv[i] << std::endl;
    return 1;
    }
  
  if (opts.input_y7 == "" || opts.input_y14 == "") {
    std::cerr << "error: both y14 and y7 ROOT source files must be supplied" << std::endl;
    return 1;
  }
  // check to make sure the input file exists
  if (!boost::filesystem::exists(opts.input_y7)) {
    std::cerr << "y7 input file does not exist: " << opts.input_y7 << std::endl;;
    return 1;
  }
  if (!boost::filesystem::exists(opts.input_y14)) {
    std::cerr << "y14 input file does not exist: " << opts.input_y14 << std::endl;;
    return 1;
  }
  
  // build output directory if it doesn't exist, using boost::filesystem
  if (opts.out_dir.empty())
    opts.out_dir = "tmp";
  boost::filesystem::path dir(opts.out_dir.c_str());
  boost::filesystem::create_directories(dir);
  
  // load root files and histograms
  TFile input_y7(opts.input_y7.c_str(), "READ");
  TFile input_y14(opts.input_y14.c_str(), "READ");
  
  TH1D* y14_pt = (TH1D*) input_y14.Get("pt");
  TH1D* y14_pt_corr = (TH1D*) input_y14.Get("ptcorr");
  TH1D* y14_refmult = (TH1D*) input_y14.Get("refmult");
  TH1D* y14_nprim = (TH1D*) input_y14.Get("nprim");
  TH1D* y14_nsel = (TH1D*) input_y14.Get("nsel");
  TH1D* y14_nhitsfit = (TH1D*) input_y14.Get("nhitsfit");
  TH1D* y14_nhitspos = (TH1D*) input_y14.Get("nhitspos");
  TH1D* y14_nhitsfitfrac = (TH1D*) input_y14.Get("nhitsfitfrac");
  TH1D* y14_dcapt = (TH1D*) input_y14.Get("dcapt");
  TProfile* y14_eff = (TProfile*) input_y14.Get("eff");
  TH2D* y14_eta_phi_0 = (TH2D*) input_y14.Get("etaphi0");
  TH2D* y14_eta_phi_1 = (TH2D*) input_y14.Get("etaphi1");
  TH2D* y14_eta_phi_2 = (TH2D*) input_y14.Get("etaphi2");
  TH2D* y14_eta_phi_3 = (TH2D*) input_y14.Get("etaphi3");

  TH1D* y7_pt = (TH1D*) input_y7.Get("pt");
  TH1D* y7_pt_corr = (TH1D*) input_y7.Get("ptcorr");
  TH1D* y7_refmult = (TH1D*) input_y7.Get("refmult");
  TH1D* y7_nprim = (TH1D*) input_y7.Get("nprim");
  TH1D* y7_nsel = (TH1D*) input_y7.Get("nsel");
  TH1D* y7_nhitsfit = (TH1D*) input_y7.Get("nhitsfit");
  TH1D* y7_nhitspos = (TH1D*) input_y7.Get("nhitspos");
  TH1D* y7_nhitsfitfrac = (TH1D*) input_y7.Get("nhitsfitfrac");
  TH1D* y7_dcapt = (TH1D*) input_y7.Get("dcapt");
  TProfile* y7_eff = (TProfile*) input_y7.Get("eff");
  TH2D* y7_eta_phi_0 = (TH2D*) input_y7.Get("etaphi0");
  TH2D* y7_eta_phi_1 = (TH2D*) input_y7.Get("etaphi1");
  TH2D* y7_eta_phi_2 = (TH2D*) input_y7.Get("etaphi2");
  TH2D* y7_eta_phi_3 = (TH2D*) input_y7.Get("etaphi3");
  
  

  
  return 0;
}
