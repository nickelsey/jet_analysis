// print_full_results.cc

#include "jet_analysis/util/root_print_routines.hh"
#include "jet_analysis/util/arg_helper.hh"

#include "TROOT.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TKey.h"
#include "TFile.h"

#include <unordered_map>
#include <iostream>

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
};

int main(int argc, char* argv[]) {
  
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
    std::cout << "AuAu input file: " << opts.auau_file;
    std::cout << " doesn't exist: exiting" << std::endl;
    return 1;
  }
  if (!boost::filesystem::exists(opts.pp_file)) {
    std::cout << "PP input file: " << opts.pp_file;
    std::cout << " doesn't exist: exiting" << std::endl;
    return 1;
  }
  
  // read in the two files
  TFile auau_file(opts.auau_file.c_str(), "READ");
  TFile pp_file(opts.pp_file.c_str(), "READ");
  
  // now we'll get the trees from the files, ignoring any objects
  // in the file that don't conform to the naming conventions from
  // the DijetWorker
  TIter auau_file_iter(auau_file.GetListOfKeys());
  TIter pp_file_iter(pp_file.GetListOfKeys());
  TKey *key;
  
  std::unordered_map<std::string, TTree*> auau_trees;
  std::unordered_map<std::string, TTree*> pp_trees;
  
  while ((key = (TKey*) auau_file_iter())) {
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
    
    auau_trees.insert({tmp, (TTree*) key->ReadObj()});
  }
  while ((key = (TKey*) pp_file_iter())) {
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
    
    pp_trees.insert({tmp, (TTree*) key->ReadObj()});
  }
  
  return 0;
}
