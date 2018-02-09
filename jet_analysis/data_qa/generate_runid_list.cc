// generate_runid_list.cc
// used to create a tree of all unique runIDs in a
// data file. Used for data analysis when plotting values
// as a function of time (runID)

#include "jet_analysis/util/arg_helper.hh"
#include "jet_analysis/util/reader_util.hh"
#include "jet_analysis/util/string_util.hh"

#include "TFile.h"
#include "TTree.h"

#include <string>
#include <set>

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

using std::string;
struct Options {
  string name        = "job";    /* output file name */
  string id          = "0";      /* job id */
  string input       = "";       /* root file/root file list*/
  string out_dir     = "";       /* directory to save output in */
  string tow_list    = "";       /* list of hot towers to remove (necessary for reader to work,
                                    even if it doesn't matter for creating a runID list) */
};

int main(int argc, char* argv[]) {
  
  // parse command line options
  Options opts;
  for (int i = 1; i < argc; ++i) {
    if (ParseStrFlag(string(argv[i]), "--name", &opts.name) ||
        ParseStrFlag(string(argv[i]), "--id", &opts.id) ||
        ParseStrFlag(string(argv[i]), "--input", &opts.input) ||
        ParseStrFlag(string(argv[i]), "--outDir", &opts.out_dir) ||
        ParseStrFlag(string(argv[i]), "--towList", &opts.tow_list)) continue;
    std::cerr << "Unknown command line option: " << argv[i] << std::endl;
    return 1;
  }

  // initialization
  // --------------
  
  // check to make sure the input file exists
  if (!boost::filesystem::exists(opts.input)) {
    std::cerr << "input file does not exist: " << opts.input << std::endl;;
    return 1;
  }
  
  // build our input chain
  TChain* chain = NewChainFromInput(opts.input);
  
  // initialize the reader(s)
  TStarJetPicoReader* reader = new TStarJetPicoReader();
  InitReaderWithDefaults(reader, chain, opts.tow_list, "");
  
  // build output directory if it doesn't exist, using boost::filesystem
  if (opts.out_dir.empty())
    opts.out_dir = "tmp";
  boost::filesystem::path dir(opts.out_dir.c_str());
  boost::filesystem::create_directories(dir);
  
  // create output file from the given directory, name & id
  string outfile_name = opts.out_dir + "/" + opts.name + opts.id + ".root";
  TFile out(outfile_name.c_str(), "RECREATE");
  
  // create the output TTree
  unsigned runID;
  TTree* runID_tree = new TTree("runid", "runid");
  runID_tree->Branch("runid", &runID);
  
  // when running we will use a std::set to hold the
  // incoming runIDs, to avoid redundancy
  std::set<unsigned> runID_set;
  
  // loop over all events
  while(reader->NextEvent()) {
    // Print out reader status every 10 seconds
    reader->PrintStatus(10);
    
    // get the run ID & the map to the index we'll use for histograms
    unsigned runID_tmp = reader->GetEvent()->GetHeader()->GetRunId();
    runID_set.insert(runID_tmp);
    
  }
  
  // now write the runIDs to the TTree, then write the TTree to file
  for (auto idx : runID_set) {
    runID = idx;
    runID_tree->Fill();
  }
  
  runID_tree->Write();
  out.Close();
  
  return 0;
}
