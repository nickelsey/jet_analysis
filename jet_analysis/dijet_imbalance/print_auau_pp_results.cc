// print_full_results.cc

#include "util/root_print_routines.hh"
#include "util/arg_helper.hh"

#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

using std::string;

// command line option results
struct Options {
  string auau_file  = "";
  string pp_file    = "";
  string out_loc    = "tmp";
  string out_prefix = "";
  int cent_low      = 0;
  int cent_high     = 8;
}

int main(int argc, char* argv[]) {
  
  Options opts;
  for (int i = 1; i < argc; ++i) {
    if (ParseStrFlag(string(argv[i]), "--auau", &opts.name) ||
        ParseStrFlag(string(argv[i]), "--pp", &opts.sj_pt) ||
        ParseStrFlag(string(argv[i]), "--outputDir", &opts.out_loc) ||
        ParseStrFlag(string(argv[i]), "--filePrefix", &opts.out_prefix) ||
        ParseIntFlag(string(argv[i]), "--centLow", &opts.cent_low) ||
        ParseIntFlag(string(argv[i]), "--centHigh", &opts.cent_high)) continue;
    std::cerr << "Unknown command line option: " << argv[i] << std::endl;
    return 1;
  }
  
  
  
  return 0;
}
