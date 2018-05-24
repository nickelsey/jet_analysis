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
  
  // load root files
  TFile input_y7(opts.input_y7.c_str(), "READ");
  TFile input_y14(opts.input_y14.c_str(), "READ");
  
//  TH2D* y14_eta_phi_0 = (TH2D*)
//  TH2D* y14_eta_phi_1 = (TH2D*)
//  TH2D* y14_eta_phi_2 = (TH2D*)
//  
//  TH2D* y7_eta_phi_0 = (TH2D*)
//  TH2D* y7_eta_phi_1 = (TH2D*)
//  TH2D* y7_eta_phi_2 = (TH2D*)
//  
  
  return 0;
}
