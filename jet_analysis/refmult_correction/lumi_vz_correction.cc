/* lumi_vz_correction.cc
 * luminosity & Vz corrections for refmult
 * distribution before glauber reweighting is done.
 */

#include "jet_analysis/util/arg_helper.hh"
#include "jet_analysis/util/string_util.hh"

#include <string>

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TF1.h"
#include "TTree.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TStyle.h"

using std::string;
struct Options {
  string input    = "";      /* input root file with refmult tree */
  string out_dir  = "tmp";   /* location to place output */
  string name     = "job";   /* name of job */
  int ref_low     = 0;       /* minimum refmult to keep when creating corrections */
  int ref_high    = 800;     /* maximum refmult to keep when creating corrections */
  double vz_low   = -30;    /* minimum event V_z */
  double vz_high  = 30;     /* maximum event V_z */
  double zdc_low  = 0;       /* minimum luminosity (ZDC coincidence rate, in Hz */
  double zdc_high = 1e5;     /* maximum luminosity (ZDC coincidence rate, in Hz */
};

int main(int argc, char* argv[]) {
  
  // Histograms will calculate gaussian errors
  // -----------------------------------------
  TH1::SetDefaultSumw2( );
  TH2::SetDefaultSumw2( );
  TH3::SetDefaultSumw2( );
  
  // set drawing preferences for histograms and graphs
  gStyle->SetOptStat(false);
  gStyle->SetOptFit(false);
  gStyle->SetOptTitle(1);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetHatchesSpacing(1.0);
  gStyle->SetHatchesLineWidth(2);
  
  Options opts;
  for (int i = 1; i < argc; ++i) {
    if (ParseStrFlag(string(argv[i]), "--input", &opts.input) ||
        ParseStrFlag(string(argv[i]), "--name", &opts.name) ||
        ParseStrFlag(string(argv[i]), "--outDir", &opts.out_dir) ||
        ParseIntFlag(string(argv[i]), "--refmultLow", &opts.ref_low) ||
        ParseIntFlag(string(argv[i]), "--refmultHigh", &opts.ref_high) ||
        ParseFloatFlag(string(argv[i]), "--VzLow", &opts.vz_low) ||
        ParseFloatFlag(string(argv[i]), "--VzHigh", &opts.vz_high) ||
        ParseFloatFlag(string(argv[i]), "--lumiLow", &opts.zdc_low) ||
        ParseFloatFlag(string(argv[i]), "--lumiHigh", &opts.zdc_high)) continue;
    std::cerr << "unknown command line option: " << argv[i] << std::endl;
    return 1;
  }
  
  // check to make sure the input file paths are sane
  if (!boost::filesystem::exists(opts.input)) {
    std::cerr << "input root file does not exist: " << opts.input << std::endl;;
    return 1;
  }
  
  // read input file
  TFile input(opts.input.c_str(), "READ");
  
  // create TTree & reader for the refmult tree
  TTree* ref_tree = (TTree*) input.Get("refMultTree");
  TTreeReader reader(ref_tree);
  TTreeReaderValue<unsigned> refmult(reader, "refMult");
  TTreeReaderValue<double> vz(reader, "vz");
  TTreeReaderValue<double> lumi(reader, "lumi");
  
  // build output directory if it doesn't exist, using boost::filesystem
  if (opts.out_dir.empty())
    opts.out_dir = "tmp";
  boost::filesystem::path dir(opts.out_dir.c_str());
  boost::filesystem::create_directories(dir);
  
  TFile out(MakeString(opts.out_dir, "/", opts.name,
                       ".root").c_str(), "RECREATE");
  
  /* first correction is done to correct the luminosity dependence
   * of the refmult distribution (lower refmult at higher luminosity)
   * and it is corrected by using a linear function, since the dependence
   * on luminosity is seen experimentally to be linear.
   */
  
  // we first fill the refmult profile as a function of luminosity
  TProfile* uncorr_lumi_prof = new TProfile("uncorr_lumi_prof", ";ZDC Rate [Hz];<refmult>",
                                            100, 0, 1e5);
  TH1D* uncorr_refmult = new TH1D("uncorr_refmult", ";refMult;counts", 800, 0, 800);
  
  while (reader.Next()) {
    if (*refmult < opts.ref_low ||
        *refmult > opts.ref_high ||
        *vz < opts.vz_low ||
        *vz > opts.vz_high ||
        *lumi < opts.zdc_low ||
        *lumi > opts.zdc_high)
      continue;
    uncorr_lumi_prof->Fill(*lumi, *refmult);
    uncorr_refmult->Fill(*refmult);
  }
  
  // fit the distribution with a straight line
  TF1* uncorr_lumi_fit = new TF1("uncorr_lumi_fit", "[0] + [1]*x", 0, 1e5);
  uncorr_lumi_fit->SetParameter(0, uncorr_lumi_prof->GetMean(2));
  uncorr_lumi_fit->SetParameter(1, 0);
  uncorr_lumi_prof->Fit(uncorr_lumi_fit);
  
  // now we will see if the correction flattens <refmult> as a function
  // of luminosity, and check that the refmult distribution is still sane.
  // if it is, we'll move onto vz corrections, so we'll fill that histogram
  // as well
  TProfile* lumi_corr_lumi_prof = new TProfile("lumi_corr_lumi_prof",
                                               ";ZDC Rate [Hz];<refmult>",
                                                100, 0, 1e5);
  TH1D* lumi_corr_refmult = new TH1D("lumi_corr_refmult", ";refMult;counts",
                                     800, 0, 800);
  
  // for vz corrections, we need vz/refmult distributions
  TH2D* lumi_corr_vz_refmult = new TH2D("lumi_corr_vz_refmult",
                                        ";V_{z} [cm];refMult",
                                        20, opts.vz_low, opts.vz_high,
                                        800, 0, 800);
  
  
  reader.Restart();
  while (reader.Next()) {
    if (*refmult < opts.ref_low ||
        *refmult > opts.ref_high ||
        *vz < opts.vz_low ||
        *vz > opts.vz_high ||
        *lumi < opts.zdc_low ||
        *lumi > opts.zdc_high)
      continue;
    double corr_refmult = *refmult * uncorr_lumi_fit->Eval(0) /
                           uncorr_lumi_fit->Eval(*lumi);
    
    lumi_corr_lumi_prof->Fill(*lumi, corr_refmult);
    lumi_corr_refmult->Fill(corr_refmult);
    lumi_corr_vz_refmult->Fill(*vz, corr_refmult);
  }
  
  // fit the corrected distribution with a straight line
  TF1* lumi_corr_lumi_fit = new TF1("lumi_corr_lumi_fit", "[0] + [1]*x", 0, 1e5);
  lumi_corr_lumi_fit->SetParameter(0, lumi_corr_lumi_prof->GetMean(2));
  lumi_corr_lumi_fit->SetParameter(1, 0);
  lumi_corr_lumi_prof->Fit(lumi_corr_lumi_fit);
  
  // make sure the correction worked before we move on
  if (lumi_corr_lumi_fit->GetParameter(1) > 1e-5) {
    std::cerr << "error: luminosity correction hasn't given zero slope, exiting" << std::endl;
    return 1;
  }
  
  // now start to correct Vz - we'll break it up into vz bins
  std::vector<TH1D*> vz_binned_refmult;
  std::vector<TF1*>  vz_binned_refmult_fit;
  
  // and we will plot the fitting parameter h as a function of V_{z}
  TH1D* vz_fit_param = new TH1D("vz_fit_param", ";V_{z};h", 20, opts.vz_low, opts.vz_high);
  for (int i = 1; i <= lumi_corr_vz_refmult->GetNbinsX(); ++i) {
    int bin = i-1;
    vz_binned_refmult.push_back((TH1D*) lumi_corr_vz_refmult
                                ->ProjectionY(MakeString("vz_bin_refmult_", i).c_str(),
                                              i, i));
    vz_binned_refmult[bin]->Scale(1.0/vz_binned_refmult[i-1]->Integral());
    
    vz_binned_refmult_fit.push_back(new TF1(MakeString("vz_bin_refmult_fit_", i).c_str(),
                                            "[0] + [0]*TMath::Erf([1]*(x-[2]))", 0, 800));
    vz_binned_refmult_fit[bin]->SetParameter(0, 0.001);
    vz_binned_refmult_fit[bin]->SetParameter(1, 0.01);
    vz_binned_refmult_fit[bin]->SetParameter(2, 565);
    
    vz_binned_refmult[bin]->Fit(vz_binned_refmult_fit[bin], "ME", "", 450, 800);
    vz_fit_param->SetBinContent(i, vz_binned_refmult_fit[bin]->GetParameter(2));
    vz_fit_param->SetBinError(i, vz_binned_refmult_fit[bin]->GetParError(2));
  }
  
  // and fit vz_fit_param with a fifth order polynomial
  TF1* vz_fit_param_fit = new TF1("vz_fit_param_fit", "pol5(0)", opts.vz_low, opts.vz_high);
  vz_fit_param_fit->SetParameters(520, 0.1, 0.1, 0.1, 0.1, 0.1);
  vz_fit_param->Fit(vz_fit_param_fit, "ME", "", opts.vz_low, opts.vz_high);
  
  
  // now we'll plot luminosity/vz/refmult and take the relevant projections
  // to check if the corrections are working
  TH3D* lumi_vz_refmult = new TH3D("lumi_vz_refmult", ";ZDC Rate[Hz];V_{z} [cm];refMult",
                                   100, 0, 1e5,
                                   20, opts.vz_low, opts.vz_high,
                                   800, 0, 800);
  reader.Restart();
  while (reader.Next()) {
    if (*refmult < opts.ref_low ||
        *refmult > opts.ref_high ||
        *vz < opts.vz_low ||
        *vz > opts.vz_high ||
        *lumi < opts.zdc_low ||
        *lumi > opts.zdc_high)
      continue;
    double corr_refmult = *refmult * uncorr_lumi_fit->Eval(0) /
                          uncorr_lumi_fit->Eval(*lumi) *
                          vz_fit_param_fit->GetParameter(0) /
                          vz_fit_param_fit->Eval(*vz);
    lumi_vz_refmult->Fill(*lumi, *vz, corr_refmult);
  }
  
  TH1D* vz_corrected_refmult = (TH1D*) lumi_vz_refmult->ProjectionZ("vz_corrected_refmult");
  TProfile* vz_corrected_lumi_prof = (TProfile*) ((TH2D*) lumi_vz_refmult->Project3D("ZX"))
                                                  ->ProfileX();
  TH2D* vz_corrected_vz_refmult = (TH2D*) lumi_vz_refmult->Project3D("ZY");
  
  std::vector<TH1D*> vz_corr_binned_refmult;
  std::vector<TF1*>  vz_corr_binned_refmult_fit;
  TH1D* vz_corr_fit_param = new TH1D("vz_corr_fit_param", ";V_{z};h", 20, opts.vz_low, opts.vz_high);
  for (int i = 1; i <= lumi_corr_vz_refmult->GetNbinsX(); ++i) {
    int bin = i-1;
    vz_corr_binned_refmult.push_back((TH1D*) vz_corrected_vz_refmult
                                ->ProjectionY(MakeString("vz_corr_bin_refmult_", i).c_str(),
                                              i, i));
    vz_corr_binned_refmult[bin]->Scale(1.0/vz_corr_binned_refmult[i-1]->Integral());
    
    vz_corr_binned_refmult_fit.push_back(new TF1(MakeString("vz_corr_bin_refmult_fit_", i).c_str(),
                                            "[0] + [0]*TMath::Erf([1]*(x-[2]))", 0, 800));
    vz_corr_binned_refmult_fit[bin]->SetParameter(0, 0.001);
    vz_corr_binned_refmult_fit[bin]->SetParameter(1, 0.01);
    vz_corr_binned_refmult_fit[bin]->SetParameter(2, 565);
    
    vz_corr_binned_refmult[bin]->Fit(vz_corr_binned_refmult_fit[bin], "ME", "", 450, 800);
    vz_corr_fit_param->SetBinContent(i, vz_corr_binned_refmult_fit[bin]->GetParameter(2));
    vz_corr_fit_param->SetBinError(i, vz_corr_binned_refmult_fit[bin]->GetParError(2));
  }
  
  // and fit vz_corr_fit_param with a fifth order polynomial
  TF1* vz_corr_fit_param_fit = new TF1("vz_corr_fit_param_fit", "[0] + [1]*x", opts.vz_low, opts.vz_high);
  vz_corr_fit_param_fit->SetParameter(0, 520);
  vz_corr_fit_param_fit->FixParameter(1, 0);
  vz_corr_fit_param->Fit(vz_corr_fit_param_fit, "ME", "", opts.vz_low, opts.vz_high);
  
    
  out.Write();
  out.Close();
  
  return 0;
}
