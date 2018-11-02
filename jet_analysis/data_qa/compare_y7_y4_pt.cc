// comapre_y7_y4_pt.cc

#include "jet_analysis/util/common.hh"
#include "jet_analysis/util/root_print_routines.hh"
#include "jet_analysis/util/string_util.hh"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TStyle.h"

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

DEFINE_string(inputRun4, "y7_y4_compare/run4_pt.root", "input run4 root file");
DEFINE_string(inputRun7, "y7_y4_compare/run7_pt.root", "input run7 root file");
DEFINE_string(outdir, "y7_y4_compare", "output directory");

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
  
  TFile y4_file(FLAGS_inputRun4.c_str(), "READ");
  TFile y7_file(FLAGS_inputRun7.c_str(), "READ");
  
  TH2D* y4_pt = (TH2D*)y4_file.Get("ptcorr");
  TH1D* y4_cent = (TH1D*)y4_file.Get("centrality");
  TH2D* y7_pt = (TH2D*)y7_file.Get("ptcorr");
  TH1D* y7_cent = ((TH2D*)y7_file.Get("refmult"))->ProjectionY("y7cent");
  
  TH1D* y4_pt_cent0 = y4_pt->ProjectionY("y4pt", 1, 1);
  y4_pt_cent0->Scale(1.0 / y4_cent->GetBinContent(1));
  TH1D* y7_pt_cent0 = y7_pt->ProjectionX("y7pt", 1, 1);
  y7_pt_cent0->RebinX(2);
  y7_pt_cent0->Scale(1.0 / y7_cent->GetBinContent(1));

  TCanvas c;
  y7_cent->Draw();
  c.SaveAs("test.pdf");

  // create our histogram and canvas options
  histogramOpts hopts;
  canvasOpts copts;
  canvasOpts coptslogz;
  coptslogz.log_z = true;
  canvasOpts coptslogy;
  coptslogy.log_y = true;
  canvasOpts cOptsBottomLeg;
  cOptsBottomLeg.leg_upper_bound = 0.4;
  cOptsBottomLeg.leg_lower_bound = 0.18;
  cOptsBottomLeg.leg_right_bound = 0.9;
  cOptsBottomLeg.leg_left_bound = 0.7;
  canvasOpts cOptsBottomLeftLeg;
  cOptsBottomLeftLeg.leg_upper_bound = 0.4;
  cOptsBottomLeftLeg.leg_lower_bound = 0.18;
  cOptsBottomLeftLeg.leg_right_bound = 0.18;
  cOptsBottomLeftLeg.leg_left_bound = 0.4;
  canvasOpts cOptsBottomLeftLegLogy;
  cOptsBottomLeftLegLogy.log_y = true;
  cOptsBottomLeftLegLogy.leg_upper_bound = 0.4;
  cOptsBottomLeftLegLogy.leg_lower_bound = 0.18;
  cOptsBottomLeftLegLogy.leg_right_bound = 0.18;
  cOptsBottomLeftLegLogy.leg_left_bound = 0.4;
  canvasOpts cOptsTopLeftLeg;
  cOptsTopLeftLeg.leg_right_bound = 0.18;
  cOptsTopLeftLeg.leg_left_bound = 0.4;

  PrintWithRatio(y4_pt_cent0, y7_pt_cent0, "y4", "y7", hopts, coptslogy, FLAGS_outdir,
                 "pt_ratio", "", "p_{T}", "1/N_{events}");

  return 0;
}