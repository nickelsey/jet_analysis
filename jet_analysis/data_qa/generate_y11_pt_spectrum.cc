// generate_y11_pt_spectrum.cc

#include "jet_analysis/util/common.hh"
#include "jet_analysis/util/trigger_lookup.hh"
#include "jet_analysis/util/reader_util.hh"
#include "jet_analysis/util/string_util.hh"
#include "jet_analysis/efficiency/run14_eff.hh"
#include "jet_analysis/efficiency/run7_eff.hh"
#include "jet_analysis/centrality/centrality_run14.hh"
#include "jet_analysis/util/root_print_routines.hh"

#include <string>
#include <iostream>
#include <cassert>
#include <vector>
#include <cmath>

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TStyle.h"

DEFINE_string(name, "y11_pt", "output file name");
DEFINE_int32(id, 0, "job id");
DEFINE_string(input, "y7_y4_compare/run4_pt.root", "input run4 root file");
DEFINE_string(outdir, "y11_pt", "output directory");

int main(int argc, char* argv[]) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);

  // set drawing preferences for histograms and graphs
  gStyle->SetOptStat(false);
  gStyle->SetOptFit(false);
  gStyle->SetLegendBorderSize(0);
  
  // build output directory if it doesn't exist, using boost::filesystem
  if (FLAGS_outdir.empty()) FLAGS_outdir = "tmp";
  boost::filesystem::path dir(FLAGS_outdir.c_str());
  boost::filesystem::create_directories(dir);
  
  // build our input chain
  TChain* chain = NewChainFromInput(FLAGS_input);
  
  // create output file from the given directory, name & id
  string outfile_name = FLAGS_outdir + "/" + FLAGS_name + std::to_string(FLAGS_id) + ".root";
  TFile out(outfile_name.c_str(), "RECREATE");
  
  // initialize the reader(s)
  TStarJetPicoReader* reader = new TStarJetPicoReader();
  InitReaderWithDefaults(reader, chain, "submit/y7_y6_bad_tower.txt", "");
  reader->GetTrackCuts()->SetDCACut(3.0);                // distance of closest approach to primary vtx
  reader->GetTrackCuts()->SetMinNFitPointsCut(20);       // minimum fit points in track reco
  reader->GetTrackCuts()->SetFitOverMaxPointsCut(0.52);  // minimum ratio of fit points used over possible
  reader->GetTrackCuts()->SetMaxPtCut(1000);             // essentially infinity - cut in eventcuts
  
  

  return 0;
}