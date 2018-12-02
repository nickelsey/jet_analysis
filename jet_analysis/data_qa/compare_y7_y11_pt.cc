// comapre_y7_y11_pt.cc

#include "jet_analysis/util/common.hh"
#include "jet_analysis/util/root_print_routines.hh"
#include "jet_analysis/util/string_util.hh"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TProfile.h"
#include "TStyle.h"

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

DEFINE_string(inputRun11, "y7_y11_compare/run11_pt.root", "input run4 root file");
DEFINE_string(inputRun7, "y7_y11_compare/run7_pt.root", "input run7 root file");
DEFINE_string(outdir, "y7_y11_compare", "output directory");

std::vector<TH1D*> ProjectXByBin(TH2D* h, string prefix, bool norm = true) {
  int bins = h->GetYaxis()->GetNbins();
  std::vector<TH1D*> ret;
  for(int i = 1; i <= bins; ++i) {
    std::string name = prefix + MakeString(i);
    ret.push_back(h->ProjectionX(name.c_str(), i, i));
    if (norm)
      ret[i-1]->Scale(1.0 / ret[i-1]->Integral());
  }
  return ret;
}

void PrintForCentralityRatio(std::vector<TH1D*> y7, std::vector<TH1D*> y11, histogramOpts hopts, canvasOpts copts,
                        string out_dir, string name, string canvas_name, string x_axis, string y_axis) {
  for (int i = 0; i < y7.size(); ++i) {
    std::string cent_name = name + MakeString(i).c_str();
    PrintWithRatio(y7[i], y11[i], "run 7", "run 14", hopts, copts, out_dir, cent_name,
                   canvas_name, x_axis, y_axis);
  }
}

template<class H>
void PrintForCentrality(std::vector<H*> y7, std::vector<H*> y11, histogramOpts hopts, canvasOpts copts,
               string out_dir, string name, string canvas_name, string x_axis, string y_axis) {
  for (int i = 0; i < y7.size(); ++i) {
    std::string cent_name = name + MakeString(i).c_str();
    Overlay1D(y7[i], y11[i], "run 7", "run 14", hopts, copts, out_dir, cent_name,
                   canvas_name, x_axis, y_axis);
  }
}


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
  
//  TH2D* y11_pt = (TH2D*)y11_file.Get("ptcorr");
//  TH1D* y11_cent = (TH1D*) ((TH2D*)y11_file.Get("refmult"))->ProjectionY("y11cent");
//  TH2D* y7_pt = (TH2D*)y7_file.Get("ptcorr");
//  TH1D* y7_cent = ((TH2D*)y7_file.Get("refmult"))->ProjectionY("y7cent");
//
//  TH1D* y11_pt_cent0 = y11_pt->ProjectionX("y11pt", 1, 1);
//  y11_pt_cent0->RebinX(2);
//  y11_pt_cent0->Scale(1.0 / y11_cent->GetBinContent(1));
//  TH1D* y7_pt_cent0 = y7_pt->ProjectionX("y7pt", 1, 1);
//  y7_pt_cent0->RebinX(2);
//  y7_pt_cent0->Scale(1.0 / y7_cent->GetBinContent(1));

    
    
//    PrintWithRatio(y11_pt_cent0, y7_pt_cent0, "y11", "y7", hopts, coptslogy, FLAGS_outdir,
//                   "pt_ratio", "", "p_{T}", "1/N_{events}");

  // create our histogram and canvas options
  histogramOpts hopts;
  canvasOpts copts;
  canvasOpts coptslogz;
  coptslogz.log_z = true;
  canvasOpts coptslogy;
  coptslogy.log_y = true;
  canvasOpts coptsBottomLeg;
  coptsBottomLeg.leg_upper_bound = 0.4;
  coptsBottomLeg.leg_lower_bound = 0.18;
  coptsBottomLeg.leg_right_bound = 0.9;
  coptsBottomLeg.leg_left_bound = 0.7;
  canvasOpts coptsBottomLeftLeg;
  coptsBottomLeftLeg.leg_upper_bound = 0.4;
  coptsBottomLeftLeg.leg_lower_bound = 0.18;
  coptsBottomLeftLeg.leg_right_bound = 0.18;
  coptsBottomLeftLeg.leg_left_bound = 0.4;
  canvasOpts coptsBottomLeftLegLogy;
  coptsBottomLeftLegLogy.log_y = true;
  coptsBottomLeftLegLogy.leg_upper_bound = 0.4;
  coptsBottomLeftLegLogy.leg_lower_bound = 0.18;
  coptsBottomLeftLegLogy.leg_right_bound = 0.18;
  coptsBottomLeftLegLogy.leg_left_bound = 0.4;
  canvasOpts coptsTopLeftLeg;
  coptsTopLeftLeg.leg_right_bound = 0.18;
  coptsTopLeftLeg.leg_left_bound = 0.4;

  TFile y11_file(FLAGS_inputRun11.c_str(), "READ");
  TFile y7_file(FLAGS_inputRun7.c_str(), "READ");

  TH2D* y11_pt = (TH2D*) y11_file.Get("pt");
  TH2D* y11_pt_corr = (TH2D*) y11_file.Get("ptcorr");
  TH2D* y11_refmult = (TH2D*) y11_file.Get("refmult");
  TH2D* y11_nprim = (TH2D*) y11_file.Get("nprim");
  TH2D* y11_nsel = (TH2D*) y11_file.Get("nsel");
  TH2D* y11_nhitsfit = (TH2D*) y11_file.Get("nhitsfit");
  TH2D* y11_nhitspos = (TH2D*) y11_file.Get("nhitspos");
  TH2D* y11_nhitsfitfrac = (TH2D*) y11_file.Get("nhitsfitfrac");
  TH3D* y11_dcapt = (TH3D*) y11_file.Get("dcapt");
  
  TH2D* y7_pt = (TH2D*) y7_file.Get("pt");
  TH2D* y7_pt_corr = (TH2D*) y7_file.Get("ptcorr");
  TH2D* y7_refmult = (TH2D*) y7_file.Get("refmult");
  TH2D* y7_nprim = (TH2D*) y7_file.Get("nprim");
  TH2D* y7_nsel = (TH2D*) y7_file.Get("nsel");
  TH2D* y7_nhitsfit = (TH2D*) y7_file.Get("nhitsfit");
  TH2D* y7_nhitspos = (TH2D*) y7_file.Get("nhitspos");
  TH2D* y7_nhitsfitfrac = (TH2D*) y7_file.Get("nhitsfitfrac");
  TH3D* y7_dcapt = (TH3D*) y7_file.Get("dcapt");
  
  std::vector<TProfile*> y11_eff(8);
  std::vector<TProfile*> y7_eff(8);
  for (int i = 0; i < 8; ++i) {
    std::string name = "eff" + MakeString(i);
    y11_eff[i] = (TProfile*) y11_file.Get(name.c_str());
    y7_eff[i] = (TProfile*) y7_file.Get(name.c_str());
  }

  // print efficiencies
  PrintForCentrality(y7_eff, y11_eff, hopts, copts, FLAGS_outdir, "efficiency", "",
                     "p_{T}", "<efficiency>");
  
  // do refmult first
  auto y7_ref_cent = ProjectXByBin(y7_refmult, "y7refmult", true);
  auto y11_ref_cent = ProjectXByBin(y11_refmult, "y11reftmult", true);
  PrintForCentrality(y7_ref_cent, y11_ref_cent, hopts, coptslogy, FLAGS_outdir, "refmult", "",
                     "refmult", "fraction");
  
  // do uncorrected pt
  
  auto y7_pt_cent = ProjectXByBin(y7_pt, "y7pt", false);
  auto y11_pt_cent = ProjectXByBin(y11_pt, "y11pt", false);
  for (int i = 0; i < y7_ref_cent.size(); ++i) {
    y7_pt_cent[i]->Scale(1.0 / y7_ref_cent[i]->GetEntries());
    y11_pt_cent[i]->Scale(1.0 / y11_ref_cent[i]->GetEntries());
  }
  PrintForCentralityRatio(y7_pt_cent, y11_pt_cent, hopts, coptslogy, FLAGS_outdir, "pt", "",
                     "p_{T}", "1/Nevents dN/dp_{T}");
  
  // and corrected pT
  auto y7_pt_corr_cent = ProjectXByBin(y7_pt_corr, "y7ptcorr", false);
  auto y11_pt_corr_cent = ProjectXByBin(y11_pt_corr, "y11ptcorr", false);
  for (int i = 0; i < y7_ref_cent.size(); ++i) {
    y7_pt_corr_cent[i]->Scale(1.0 / y7_ref_cent[i]->GetEntries());
    y11_pt_corr_cent[i]->Scale(1.0 / y11_ref_cent[i]->GetEntries());
    
    for (int j = 1; j <= y11_pt_corr_cent[i]->GetNbinsX(); ++j) {
      if (!std::isfinite(y11_pt_corr_cent[i]->GetBinContent(j))) {
        y11_pt_corr_cent[i]->SetBinContent(j, 0.00);
        y11_pt_corr_cent[i]->SetBinError(j, 0.00);
      }
    }
  }
  
    
  PrintForCentralityRatio(y7_pt_corr_cent, y11_pt_corr_cent, hopts, coptslogy, FLAGS_outdir, "ptcorr", "",
                          "p_{T}", "1/Nevents dN/dp_{T}");
  
  
  auto y7_nprim_cent = ProjectXByBin(y7_nprim, "y7nprim", true);
  auto y11_nprim_cent = ProjectXByBin(y11_nprim, "y11nprim", true);
  PrintForCentrality(y7_nprim_cent, y11_nprim_cent, hopts, coptslogy, FLAGS_outdir, "nprim", "",
                     "N_{primaries}", "fraction");
  
  auto y7_nhit_cent = ProjectXByBin(y7_nhitsfit, "y7nhitsfit", true);
  auto y11_nhit_cent = ProjectXByBin(y11_nhitsfit, "y11nhitsfit", true);
  PrintForCentrality(y7_nhit_cent, y11_nhit_cent, hopts, coptsBottomLeftLeg, FLAGS_outdir, "nhit", "",
                     "fit points", "fraction");
  
  auto y7_nhitpos_cent = ProjectXByBin(y7_nhitspos, "y7nhitspos", true);
  auto y11_nhitpos_cent = ProjectXByBin(y11_nhitspos, "y11nhitspos", true);
  PrintForCentrality(y7_nhitpos_cent, y11_nhitpos_cent, hopts, coptsBottomLeftLeg, FLAGS_outdir, "nhitspos", "",
                     "fit points possible", "fraction");
  
  auto y7_nhitfrac_cent = ProjectXByBin(y7_nhitsfitfrac, "y7nhitsfrac", true);
  auto y11_nhitfrac_cent = ProjectXByBin(y11_nhitsfitfrac, "y11nhitsfrac", true);
  PrintForCentrality(y7_nhitfrac_cent, y11_nhitfrac_cent, hopts, coptsBottomLeftLeg, FLAGS_outdir, "nhitsfrac", "",
                     "fit points / possible", "fraction");
  
  
  TFile out(MakeString(FLAGS_outdir, "/ratio.root").c_str(), "RECREATE");
  for (int i = 0; i < y7_pt_corr_cent.size(); ++i) {
    TProfile* y11_eff_proj = y11_eff[i];
    TProfile* y7_eff_proj = y7_eff[i];
    std::string cent_ratio_name = "eff_ratio_" + std::to_string(i);
    TH1D* ratio = (TH1D*) y7_pt_corr_cent[i]->Clone(cent_ratio_name.c_str());
    if (ratio->GetEntries() <= 0)
      continue;
    ratio->Divide(y11_pt_corr_cent[i]);
    
    TF1* fit = new TF1(MakeString("fit_", i).c_str(),"[0]", 2, 5);
    fit->SetParameter(0, 1.3);
    ratio->Fit(fit, "", "", 2.0, 5.0);
    
    out.cd();
    ratio->Write();
  }
  
  out.Close();
  
  
  return 0;
}