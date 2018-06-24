#include "jet_analysis/util/arg_helper.hh"
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

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TStyle.h"

using std::string;

struct Options {
  string input_y7    = "";       /* root file for y7 data*/
  string input_y14   = "";       /* root file for y14 data */
  string out_dir     = "tmp";    /* directory to save output in */
};


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

void PrintForCentralityRatio(std::vector<TH1D*> y7, std::vector<TH1D*> y14, histogramOpts hOpts, canvasOpts cOpts,
                        string out_dir, string name, string canvas_name, string x_axis, string y_axis) {
  for (int i = 0; i < y7.size(); ++i) {
    std::string cent_name = name + MakeString(i).c_str();
    PrintWithRatio(y7[i], y14[i], "run 7", "run 14", hOpts, cOpts, out_dir, cent_name,
                   canvas_name, x_axis, y_axis);
  }
}

void PrintForCentrality(std::vector<TH1D*> y7, std::vector<TH1D*> y14, histogramOpts hOpts, canvasOpts cOpts,
               string out_dir, string name, string canvas_name, string x_axis, string y_axis) {
  for (int i = 0; i < y7.size(); ++i) {
    std::string cent_name = name + MakeString(i).c_str();
    Overlay1D(y7[i], y14[i], "run 7", "run 14", hOpts, cOpts, out_dir, cent_name,
                   canvas_name, x_axis, y_axis);
  }
}


int main(int argc, char* argv[]) {
  
  // set drawing preferences for histograms and graphs
  gStyle->SetOptStat(false);
  gStyle->SetOptFit(false);
  gStyle->SetLegendBorderSize(0);
  
  // create histogram options
  histogramOpts hOpts;
  canvasOpts cOpts;
  canvasOpts cOptsNoLeg;
  cOptsNoLeg.do_legend = false;
  canvasOpts cOptsLogy;
  cOptsLogy.log_y = true;
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
  
  TH2D* y14_pt = (TH2D*) input_y14.Get("pt");
  TH2D* y14_pt_corr = (TH2D*) input_y14.Get("ptcorr");
  TH2D* y14_refmult = (TH2D*) input_y14.Get("refmult");
  TH2D* y14_nprim = (TH2D*) input_y14.Get("nprim");
  TH2D* y14_nsel = (TH2D*) input_y14.Get("nsel");
  TH2D* y14_nhitsfit = (TH2D*) input_y14.Get("nhitsfit");
  TH2D* y14_nhitspos = (TH2D*) input_y14.Get("nhitspos");
  TH2D* y14_nhitsfitfrac = (TH2D*) input_y14.Get("nhitsfitfrac");
  TH3D* y14_dcapt = (TH3D*) input_y14.Get("dcapt");
  TH2D* y14_eta_phi_0 = (TH2D*) input_y14.Get("etaphi0");
  TH2D* y14_eta_phi_1 = (TH2D*) input_y14.Get("etaphi1");
  TH2D* y14_eta_phi_2 = (TH2D*) input_y14.Get("etaphi2");
  TH2D* y14_eta_phi_3 = (TH2D*) input_y14.Get("etaphi3");
  
  TH2D* y7_pt = (TH2D*) input_y7.Get("pt");
  TH2D* y7_pt_corr = (TH2D*) input_y7.Get("ptcorr");
  TH2D* y7_refmult = (TH2D*) input_y7.Get("refmult");
  TH2D* y7_nprim = (TH2D*) input_y7.Get("nprim");
  TH2D* y7_nsel = (TH2D*) input_y7.Get("nsel");
  TH2D* y7_nhitsfit = (TH2D*) input_y7.Get("nhitsfit");
  TH2D* y7_nhitspos = (TH2D*) input_y7.Get("nhitspos");
  TH2D* y7_nhitsfitfrac = (TH2D*) input_y7.Get("nhitsfitfrac");
  TH3D* y7_dcapt = (TH3D*) input_y7.Get("dcapt");
  TH2D* y7_eta_phi_0 = (TH2D*) input_y7.Get("etaphi0");
  TH2D* y7_eta_phi_1 = (TH2D*) input_y7.Get("etaphi1");
  TH2D* y7_eta_phi_2 = (TH2D*) input_y7.Get("etaphi2");
  TH2D* y7_eta_phi_3 = (TH2D*) input_y7.Get("etaphi3");
  
  std::vector<TProfile*> y14_eff(8);
  std::vector<TProfile*> y7_eff(8);
  for (int i = 0; i < 8; ++i) {
    std::string name = "eff" + MakeString(i);
    y14_eff[i] = (TProfile*) input_y14.Get(name.c_str());
    y7_eff[i] = (TProfile*) input_y7.Get(name.c_str());
  }
  
  // do refmult first
  auto y7_ref_cent = ProjectXByBin(y7_refmult, "y7refmult", true);
  auto y14_ref_cent = ProjectXByBin(y14_refmult, "y14reftmult", true);
  PrintForCentrality(y7_ref_cent, y14_ref_cent, hOpts, cOptsLogy, opts.out_dir, "refmult", "",
                     "refmult", "fraction");
  
  // do uncorrected pt
  
  auto y7_pt_cent = ProjectXByBin(y7_pt, "y7pt", false);
  auto y14_pt_cent = ProjectXByBin(y14_pt, "y14pt", false);
  for (int i = 0; i < y7_ref_cent.size(); ++i) {
    y7_pt_cent[i]->Scale(1.0 / y7_ref_cent[i]->GetEntries());
    y14_pt_cent[i]->Scale(1.0 / y14_ref_cent[i]->GetEntries());
  }
  PrintForCentralityRatio(y7_pt_cent, y14_pt_cent, hOpts, cOptsLogy, opts.out_dir, "pt", "",
                     "p_{T}", "1/Nevents dN/dp_{T}");
  
  auto y7_nprim_cent = ProjectXByBin(y7_nprim, "y7nprim", true);
  auto y14_nprim_cent = ProjectXByBin(y14_nprim, "y14nprim", true);
  PrintForCentrality(y7_nprim_cent, y14_nprim_cent, hOpts, cOptsLogy, opts.out_dir, "nprim", "",
                     "N_{primaries}", "fraction");
  
  auto y7_nhit_cent = ProjectXByBin(y7_nhitsfit, "y7nhitsfit", true);
  auto y14_nhit_cent = ProjectXByBin(y14_nhitsfit, "y14nhitsfit", true);
  PrintForCentrality(y7_nhit_cent, y14_nhit_cent, hOpts, cOptsBottomLeftLeg, opts.out_dir, "nhit", "",
                     "fit points", "fraction");
  
  auto y7_nhitpos_cent = ProjectXByBin(y7_nhitspos, "y7nhitspos", true);
  auto y14_nhitpos_cent = ProjectXByBin(y14_nhitspos, "y14nhitspos", true);
  PrintForCentrality(y7_nhitpos_cent, y14_nhitpos_cent, hOpts, cOptsBottomLeftLeg, opts.out_dir, "nhitspos", "",
                     "fit points possible", "fraction");
  
  auto y7_nhitfrac_cent = ProjectXByBin(y7_nhitsfitfrac, "y7nhitsfrac", true);
  auto y14_nhitfrac_cent = ProjectXByBin(y14_nhitsfitfrac, "y14nhitsfrac", true);
  PrintForCentrality(y7_nhitfrac_cent, y14_nhitfrac_cent, hOpts, cOptsBottomLeftLeg, opts.out_dir, "nhitsfrac", "",
                     "fit points / possible", "fraction");
  
  
  
//  // normalize
//  y14_eta_phi_0->Scale(1.0 / y14_refmult->GetEntries());
//  y14_eta_phi_1->Scale(1.0 / y14_refmult->GetEntries());
//  y14_eta_phi_2->Scale(1.0 / y14_refmult->GetEntries());
//  y14_eta_phi_3->Scale(1.0 / y14_refmult->GetEntries());
//  y14_pt->Scale(1.0 / y14_refmult->GetEntries());
//  y14_pt_corr->Scale(1.0 / y14_refmult->GetEntries());
//  y14_refmult->Scale(1.0 / y14_refmult->Integral());
//  y14_nprim->Scale(1.0 / y14_nprim->Integral());
//  y14_nsel->Scale(1.0 / y14_nsel->Integral());
//  y14_nhitsfit->Scale(1.0 / y14_nhitsfit->Integral());
//  y14_nhitsfitfrac->Scale(1.0 / y14_nhitsfitfrac->Integral());
//  y14_nhitspos->Scale(1.0 / y14_nhitspos->Integral());
//  y14_dcapt->Scale(1.0 / y14_dcapt->Integral());
//
//  y7_eta_phi_0->Scale(1.0 / y7_refmult->GetEntries());
//  y7_eta_phi_1->Scale(1.0 / y7_refmult->GetEntries());
//  y7_eta_phi_2->Scale(1.0 / y7_refmult->GetEntries());
//  y7_eta_phi_3->Scale(1.0 / y7_refmult->GetEntries());
//  y7_pt->Scale(1.0 / y7_refmult->GetEntries());
//  y7_pt_corr->Scale(1.0 / y7_refmult->GetEntries());
//  y7_refmult->Scale(1.0 / y7_refmult->Integral());
//  y7_nprim->Scale(1.0 / y7_nprim->Integral());
//  y7_nsel->Scale(1.0 / y7_nsel->Integral());
//  y7_nhitsfit->Scale(1.0 / y7_nhitsfit->Integral());
//  y7_nhitsfitfrac->Scale(1.0 / y7_nhitsfitfrac->Integral());
//  y7_nhitspos->Scale(1.0 / y7_nhitspos->Integral());
//  y7_dcapt->Scale(1.0 / y7_dcapt->Integral());
//
//  TH1D* y14_dca = (TH1D*) y14_dcapt->ProjectionY();
//  y14_dca->SetName("y14dca");
//  TH1D* y7_dca = (TH1D*) y7_dcapt->ProjectionY();
//  y7_dca->SetName("y7dca");
//  y14_dca->Scale(1.0 / y14_dca->Integral());
//  y7_dca->Scale(1.0 / y7_dca->Integral());
//
//  PrintWithRatio(y7_pt, y14_pt, "run 7", "run 14", hOpts, cOptsLogy, opts.out_dir, "uncorr_pt",
//            "", "p_{T}", "1/Nevents dN/dp_{T}");
//  PrintWithRatio(y7_pt_corr, y14_pt_corr, "run 7", "run 14", hOpts, cOptsLogy, opts.out_dir, "corr_pt",
//            "", "p_{T}", "1/Nevents dN/dp_{T}");
//  Overlay1D(y7_refmult, y14_refmult, "run 7", "run 14", hOpts, cOptsBottomLeftLegLogy, opts.out_dir, "refmult",
//            "", "refmult", "fraction");
//  Overlay1D(y7_nprim, y14_nprim, "run 7", "run 14", hOpts, cOptsLogy, opts.out_dir, "prim", "",
//            "N_{primaries}", "fraction");
//  Overlay1D(y7_nsel, y14_nsel, "run 7", "run 14", hOpts, cOptsLogy, opts.out_dir, "nsel", "",
//            "N_{tracks}", "fraction");
//  Overlay1D(y7_nhitsfit, y14_nhitsfit, "run 7", "run 14", hOpts, cOptsBottomLeftLeg, opts.out_dir, "nhits",
//            "", "nhits fit", "fraction");
//  Overlay1D(y7_nhitspos, y14_nhitspos, "run 7", "run 14", hOpts, cOptsBottomLeftLeg, opts.out_dir, "nhitspos",
//            "", "nhits possible", "fraction");
//  Overlay1D(y7_nhitsfitfrac, y14_nhitsfitfrac, "run 7", "run 14", hOpts, cOptsBottomLeftLeg, opts.out_dir, "nhitsfitfrac",
//            "", "nhits fit fraction", "fraction");
//  Overlay1D(y7_dca, y14_dca, "run 7", "run 14", hOpts, cOptsLogy, opts.out_dir, "dca", "", "DCA [cm]", "fraction");
//
//  Print2DSimple(y7_eta_phi_0, hOpts, cOpts, opts.out_dir, "y7_ep_0", "", "#eta", "#phi");
//  Print2DSimple(y7_eta_phi_1, hOpts, cOpts, opts.out_dir, "y7_ep_1", "", "#eta", "#phi");
//  Print2DSimple(y7_eta_phi_2, hOpts, cOpts, opts.out_dir, "y7_ep_2", "", "#eta", "#phi");
//  Print2DSimple(y7_eta_phi_3, hOpts, cOpts, opts.out_dir, "y7_ep_3", "", "#eta", "#phi");
//
//  Print2DSimple(y14_eta_phi_0, hOpts, cOpts, opts.out_dir, "y14_ep_0", "", "#eta", "#phi");
//  Print2DSimple(y14_eta_phi_1, hOpts, cOpts, opts.out_dir, "y14_ep_1", "", "#eta", "#phi");
//  Print2DSimple(y14_eta_phi_2, hOpts, cOpts, opts.out_dir, "y14_ep_2", "", "#eta", "#phi");
//  Print2DSimple(y14_eta_phi_3, hOpts, cOpts, opts.out_dir, "y14_ep_3", "", "#eta", "#phi");
//
//  // divide
//  y14_eta_phi_0->Divide(y7_eta_phi_0);
//  y14_eta_phi_1->Divide(y7_eta_phi_1);
//  y14_eta_phi_2->Divide(y7_eta_phi_2);
//  y14_eta_phi_3->Divide(y7_eta_phi_3);
//
//  Print2DSimple(y14_eta_phi_0, hOpts, cOpts, opts.out_dir, "y14_y7_0", "", "#eta", "#phi");
//  Print2DSimple(y14_eta_phi_1, hOpts, cOpts, opts.out_dir, "y14_y7_1", "", "#eta", "#phi");
//  Print2DSimple(y14_eta_phi_2, hOpts, cOpts, opts.out_dir, "y14_y7_2", "", "#eta", "#phi");
//  Print2DSimple(y14_eta_phi_3, hOpts, cOpts, opts.out_dir, "y14_y7_3", "", "#eta", "#phi");
//
  return 0;
}

