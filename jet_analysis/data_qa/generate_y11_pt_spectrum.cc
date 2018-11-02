// generate_y11_pt_spectrum.cc

#include "jet_analysis/util/common.hh"
#include "jet_analysis/util/trigger_lookup.hh"
#include "jet_analysis/util/reader_util.hh"
#include "jet_analysis/util/string_util.hh"
#include "jet_analysis/efficiency/run11_eff.hh"
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
DEFINE_string(effFile, "submit/y11_effic.root", "efficiency curves for Run 11");
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
  
  // load Run11 efficiency curves
  Run11Eff eff(FLAGS_effFile);

  // initialize the reader(s)
  TStarJetPicoReader* reader = new TStarJetPicoReader();
  InitReaderWithDefaults(reader, chain, "submit/y7_y6_bad_tower.txt", "");
  reader->GetTrackCuts()->SetDCACut(3.0);                // distance of closest approach to primary vtx
  reader->GetTrackCuts()->SetMinNFitPointsCut(20);       // minimum fit points in track reco
  reader->GetTrackCuts()->SetFitOverMaxPointsCut(0.52);  // minimum ratio of fit points used over possible
  reader->GetTrackCuts()->SetMaxPtCut(1000);             // essentially infinity - cut in eventcuts
  
  // create output file from the given directory, name & id
  string outfile_name = FLAGS_outdir + "/" + FLAGS_name + std::to_string(FLAGS_id) + ".root";
  TFile out(outfile_name.c_str(), "RECREATE");

  // Histograms will calculate gaussian errors
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();
  
  int cent_bins = 9;

  // create pT histogram
  TH2D* nglobnprim = new TH2D("nglobnprim", "", 500, 0, 5000, 500, 0, 1500);
  TH2D* refmultnprim = new TH2D("refmultnprim", "", 800, 0, 800, 500, 0, 1500);
  TH2D* pt = new TH2D("pt", ";p_{T};centrality", 100, 0, 5, cent_bins, -0.5, cent_bins - 0.5);
  TH2D* pt_corr = new TH2D("ptcorr", ";p_{T};centrality", 100, 0, 5,cent_bins, -0.5, cent_bins - 0.5);
  TH2D* refmult = new TH2D("refmult", ";refmult;centrality", 800, 0, 800, cent_bins, -0.5, cent_bins - 0.5);
  TH2D* frac = new TH2D("discarded", "", 10, 0, 1.0, cent_bins, -0.5, cent_bins - 0.5);
  TH2D* nprim = new TH2D("nprim", "", 100, 0, 1200, cent_bins, -0.5, cent_bins - 0.5);
  TH2D* nsel = new TH2D("nsel", "", 100, 0, 12000, cent_bins, -0.5, cent_bins - 0.5);
  TH2D* nhitsfit = new TH2D("nhitsfit", "", 50, 0, 50, cent_bins, -0.5, cent_bins - 0.5);
  TH2D* nhitspos = new TH2D("nhitspos", "", 50, 0, 50, cent_bins, -0.5, cent_bins - 0.5);
  TH2D* nhitsfitfrac = new TH2D("nhitsfitfrac", "", 50, 0, 1.0, cent_bins, -0.5, cent_bins - 0.5);
  TH3D* dca = new TH3D("dcapt", "", 50, 0, 3, 50, 0, 5, cent_bins, -0.5, cent_bins - 0.5);
  
  

  return 0;
}