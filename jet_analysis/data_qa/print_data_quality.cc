// print_data_quality.cc

#include "jet_analysis/util/arg_helper.hh"
#include "jet_analysis/util/string_util.hh"
#include "jet_analysis/util/root_print_routines.hh"
#include "jet_analysis/util/histogram_routines.hh"

#include <string>
#include <unordered_map>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TF1.h"
#include "TStyle.h"

using std::string;

struct Options {
  string input       = "";       /* root file/root file list*/
  string out_dir     = "tmp";    /* directory to save output in */
  bool   use_y7      = false;    /* whether or not to include y7 comparisons */
  bool   use_high    = true;     /* flag to include high luminosity y14 data */
  bool   use_mid     = true;     /* flag to include mid luminosity y14 data */
  bool   use_low     = true;     /* flag to include low luminosity y14 data */
  bool   use_pre     = false;    /* flag to include AuAu_200_production_2014 data */
};

int main(int argc, char* argv[]) {
  
  // set drawing preferences for histograms and graphs
  gStyle->SetOptStat(false);
  gStyle->SetOptFit(false);
  gStyle->SetLegendBorderSize(0);
  
  std::vector<string> prefixes{"y7", "pre", "low", "mid", "high"};
  
  // create histogram options
  histogramOpts hOpts;
  canvasOpts cOpts;
  canvasOpts cOptsNoLeg;
  cOptsNoLeg.do_legend = false;
  
  // parse command line options
  // --------------------------
  Options opts;
  for (int i = 1; i < argc; ++i) {
    if (ParseStrFlag(string(argv[i]), "--input", &opts.input) ||
        ParseStrFlag(string(argv[i]), "--outDir", &opts.out_dir) ||
        ParseBoolFlag(string(argv[i]), "--y7", &opts.use_y7) ||
        ParseBoolFlag(string(argv[i]), "--high", &opts.use_high) ||
        ParseBoolFlag(string(argv[i]), "--mid", &opts.use_mid) ||
        ParseBoolFlag(string(argv[i]), "--low", &opts.use_low) ||
        ParseBoolFlag(string(argv[i]), "--pre", &opts.use_pre)) continue;
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
  
  // build output directory if it doesn't exist, using boost::filesystem
  if (opts.out_dir.empty())
    opts.out_dir = "tmp";
  boost::filesystem::path dir(opts.out_dir.c_str());
  boost::filesystem::create_directories(dir);
  
  // input file
  TFile input(opts.input.c_str(), "READ");
  
  // we will store histograms in a dictionary for
  // easy lookup
  std::unordered_map<string, TH2D*> runid_refmult;
  std::unordered_map<string, TH2D*> runid_grefmult;
  std::unordered_map<string, TH2D*> runid_nprim;
  std::unordered_map<string, TH2D*> runid_nglob;
  std::unordered_map<string, TH2D*> zdc_refmult;
  std::unordered_map<string, TH2D*> bbc_refmult;
  std::unordered_map<string, TH2D*> zdc_grefmult;
  std::unordered_map<string, TH2D*> vz_refmult;
  std::unordered_map<string, TH2D*> refmult_grefmult;
  std::unordered_map<string, TH2D*> prim_glob;
  std::unordered_map<string, TH2D*> runid_vx;
  std::unordered_map<string, TH2D*> runid_vy;
  std::unordered_map<string, TH2D*> runid_vz;
  std::unordered_map<string, TH2D*> vx_vy;
  std::unordered_map<string, TH2D*> vx_vz;
  std::unordered_map<string, TH2D*> vy_vz;
  std::unordered_map<string, TH2D*> zdc_vz;
  std::unordered_map<string, TH2D*> zdc_vzvpdvz;
  std::unordered_map<string, TH2D*> runid_zdc;
  std::unordered_map<string, TH2D*> runid_bbc;
  std::unordered_map<string, TH2D*> runid_px;
  std::unordered_map<string, TH2D*> runid_py;
  std::unordered_map<string, TH2D*> runid_pz;
  std::unordered_map<string, TH2D*> runid_pt;
  std::unordered_map<string, TH2D*> zdc_px;
  std::unordered_map<string, TH2D*> zdc_py;
  std::unordered_map<string, TH2D*> zdc_pz;
  std::unordered_map<string, TH2D*> zdc_pt;
  std::unordered_map<string, TH2D*> runid_dca;
  std::unordered_map<string, TH2D*> runid_fit;
  std::unordered_map<string, TH2D*> runid_eta;
  std::unordered_map<string, TH2D*> runid_phi;
  std::unordered_map<string, TH2D*> track_etaphi;
  std::unordered_map<string, TH2D*> runid_towe;
  std::unordered_map<string, TH2D*> runid_towet;
  std::unordered_map<string, TH2D*> runid_towadc;
  std::unordered_map<string, TH2D*> zdc_towe;
  std::unordered_map<string, TH2D*> zdc_towet;
  std::unordered_map<string, TH2D*> zdc_towadc;
  std::unordered_map<string, TH2D*> tow_towe;
  std::unordered_map<string, TH2D*> tow_towet;
  std::unordered_map<string, TH2D*> tow_towadc;
  std::unordered_map<string, TH2D*> tow_etaphi;
  
  // read histograms in from input
  for (auto prefix : prefixes) {
    runid_refmult.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "runidrefmult").c_str())});
    runid_grefmult.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "runidgrefmult").c_str())});
    runid_nprim.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "runidnprim").c_str())});
    runid_nglob.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "runidnglob").c_str())});
    zdc_refmult.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "zdcrefmult").c_str())});
    bbc_refmult.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "bbcrefmult").c_str())});
    zdc_grefmult.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "zdcgrefmult").c_str())});
    vz_refmult.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "vzrefmult").c_str())});
    refmult_grefmult.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "refgrefmult").c_str())});
    prim_glob.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "primglob").c_str())});
    runid_vx.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "runidvx").c_str())});
    runid_vy.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "runidvy").c_str())});
    runid_vz.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "runidvz").c_str())});
    vx_vy.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "vxvy").c_str())});
    vx_vz.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "vxvz").c_str())});
    vy_vz.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "vyvz").c_str())});
    zdc_vz.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "zdcvz").c_str())});
    zdc_vzvpdvz.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "vzvpdvz").c_str())});
    runid_zdc.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "runidzdc").c_str())});
    runid_bbc.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "runidbbc").c_str())});
    runid_px.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "runidpx").c_str())});
    runid_py.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "runidpy").c_str())});
    runid_pz.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "runidpz").c_str())});
    runid_pt.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "runidpt").c_str())});
    zdc_px.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "zdcpx").c_str())});
    zdc_py.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "zdcpy").c_str())});
    zdc_pz.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "zdcpz").c_str())});
    zdc_pt.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "zdcpt").c_str())});
    runid_dca.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "runiddca").c_str())});
    runid_fit.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "runidfit").c_str())});
    runid_eta.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "runideta").c_str())});
    runid_phi.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "runidphi").c_str())});
    track_etaphi.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "tracketaphi").c_str())});
    runid_towe.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "runidtowe").c_str())});
    runid_towet.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "runidtowet").c_str())});
    runid_towadc.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "runidtowadc").c_str())});
    zdc_towe.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "zdctowe").c_str())});
    zdc_towet.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "zdctowet").c_str())});
    zdc_towadc.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "zdctowadc").c_str())});
    tow_towe.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "towtowe").c_str())});
    tow_towet.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "towtowet").c_str())});
    tow_towadc.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "towtowadc").c_str())});
    tow_etaphi.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "towetaphi").c_str())});
  }
  
  
  // printing routines
  // -----------------
  
  // first, nprimary vs runid
  std::vector<TProfile*> runid_nprim_avg_set;
  std::vector<std::string> runid_nprim_avg_prefix;
  if (opts.use_high) {
    runid_nprim_avg_set.push_back((TProfile*) runid_nprim["high"]->ProfileX());
    runid_nprim_avg_prefix.push_back("high");
  }
  if (opts.use_low) {
    runid_nprim_avg_set.push_back((TProfile*) runid_nprim["low"]->ProfileX());
    runid_nprim_avg_prefix.push_back("low");
  }
  if (opts.use_mid) {
    runid_nprim_avg_set.push_back((TProfile*) runid_nprim["mid"]->ProfileX());
    runid_nprim_avg_prefix.push_back("mid");
  }
  if (opts.use_pre) {
    runid_nprim_avg_set.push_back((TProfile*) runid_nprim["pre"]->ProfileX());
    runid_nprim_avg_prefix.push_back("pre");
  }
  runid_nprim_avg_set[0]->GetXaxis()->SetRangeUser(1200, 1400);
  runid_nprim_avg_set[0]->GetYaxis()->SetRangeUser(200, 500);
  Overlay1D(runid_nprim_avg_set, runid_nprim_avg_prefix, hOpts, cOpts, opts.out_dir, "runidnprim",
            "", "run ID", "<N_{primary}>", "", false);
  
  // then, nglobal vs runid
  std::vector<TProfile*> runid_nglob_avg_set;
  std::vector<std::string> runid_nglob_avg_prefix;
  if (opts.use_high) {
    runid_nglob_avg_set.push_back((TProfile*) runid_nglob["high"]->ProfileX());
    runid_nglob_avg_prefix.push_back("high");
  }
  if (opts.use_low) {
    runid_nglob_avg_set.push_back((TProfile*) runid_nglob["low"]->ProfileX());
    runid_nglob_avg_prefix.push_back("low");
  }
  if (opts.use_mid) {
    runid_nglob_avg_set.push_back((TProfile*) runid_nglob["mid"]->ProfileX());
    runid_nglob_avg_prefix.push_back("mid");
  }
  if (opts.use_pre) {
    runid_nglob_avg_set.push_back((TProfile*) runid_nglob["pre"]->ProfileX());
    runid_nglob_avg_prefix.push_back("pre");
  }
  runid_nglob_avg_set[0]->GetXaxis()->SetRangeUser(1200, 1400);
  runid_nglob_avg_set[0]->GetYaxis()->SetRangeUser(1000, 2500);
  Overlay1D(runid_nglob_avg_set, runid_nglob_avg_prefix, hOpts, cOpts, opts.out_dir, "runidnglob",
            "", "run ID", "<N_{global}>", "", false);
  
//  // event level
//  // -----------
//  TH2D* lumi_refmult = (TH2D*) zdc_refmult["pre"]->Clone("zdc_refmult_clone");
//  lumi_refmult->Add(zdc_refmult["low"]);
//  lumi_refmult->Add(zdc_refmult["mid"]);
//  lumi_refmult->Add(zdc_refmult["high"]);
//  TH1D* low_lumi_refmult = (TH1D*) lumi_refmult->ProjectionY("_low", 1, 33);
//  low_lumi_refmult->Scale(1.0/low_lumi_refmult->Integral());
//  TH1D* mid_lumi_refmult = (TH1D*) lumi_refmult->ProjectionY("_mid", 34, 66);
//  mid_lumi_refmult->Scale(1.0/mid_lumi_refmult->Integral());
//  TH1D* high_lumi_refmult = (TH1D*) lumi_refmult->ProjectionY("_high", 67, 100);
//  high_lumi_refmult->Scale(1.0/high_lumi_refmult->Integral());
//  std::vector<TH1D*> lumi_refmult_bins{low_lumi_refmult, mid_lumi_refmult, high_lumi_refmult};
//  std::vector<string> lumi_refmult_bin_names{"0-33 kHz", "34-66 kHz", "67-100 kHz"};
//  Overlay1D(lumi_refmult_bins, lumi_refmult_bin_names, hOpts, cOpts, opts.out_dir, "lumirefmult", "",
//            "refMult", "fraction", "ZDC Rate");
//
//  // print refmult as a function of vz
//  TH2D* vz_ref = (TH2D*) vz_refmult["pre"]->Clone("vz_ref_clone");
//  vz_ref->Add(vz_refmult["low"]);
//  vz_ref->Add(vz_refmult["mid"]);
//  vz_ref->Add(vz_refmult["high"]);
//  TH1D* low_vz_refmult = (TH1D*) vz_ref->ProjectionY("_low", 1, 33);
//  low_vz_refmult->Scale(1.0/low_vz_refmult->Integral());
//  TH1D* mid_vz_refmult = (TH1D*) vz_ref->ProjectionY("_mid", 34, 66);
//  mid_vz_refmult->Scale(1.0/mid_vz_refmult->Integral());
//  TH1D* high_vz_refmult = (TH1D*) vz_ref->ProjectionY("_high", 67, 100);
//  high_vz_refmult->Scale(1.0/high_vz_refmult->Integral());
//  std::vector<TH1D*> vz_ref_bins{low_vz_refmult, mid_vz_refmult, high_vz_refmult};
//  std::vector<string> vz_ref_bin_names{"-30 < V_{z} < -10 cm",
//                                       "-10 < V_{z} < 10 cm",
//                                       "10 < V_{z} < 30 cm"};
//  Overlay1D(vz_ref_bins, vz_ref_bin_names, hOpts, cOpts, opts.out_dir, "vzrefmult", "",
//            "refMult", "fraction", "V_{z}");
//
//  // vz, vx, vy
//  TH2D* runid_vz_full = (TH2D*) runid_vz["pre"]->Clone("vz_full_clone");
//  runid_vz_full->Add(runid_vz["low"]);
//  runid_vz_full->Add(runid_vz["mid"]);
//  runid_vz_full->Add(runid_vz["high"]);
//  TProfile* vz_prof = (TProfile*) runid_vz_full->ProfileX("vz_runid_profile");
//  PrettyPrint1D(vz_prof, hOpts, cOptsNoLeg, "", opts.out_dir, "run_avg_vz", "", "runID", "<V_{z}>");
//
//  TH2D* runid_vx_full = (TH2D*) runid_vx["pre"]->Clone("vx_full_clone");
//  runid_vx_full->Add(runid_vx["low"]);
//  runid_vx_full->Add(runid_vx["mid"]);
//  runid_vx_full->Add(runid_vx["high"]);
//  TProfile* vx_prof = (TProfile*) runid_vx_full->ProfileX("vx_runid_profile");
//  PrettyPrint1D(vx_prof, hOpts, cOptsNoLeg, "", opts.out_dir, "run_avg_vx", "", "runID", "<V_{x}>");
//
//  TH2D* runid_vy_full = (TH2D*) runid_vy["pre"]->Clone("vy_full_clone");
//  runid_vy_full->Add(runid_vy["low"]);
//  runid_vy_full->Add(runid_vy["mid"]);
//  runid_vy_full->Add(runid_vy["high"]);
//  TProfile* vy_prof = (TProfile*) runid_vy_full->ProfileX("vy_runid_profile");
//  PrettyPrint1D(vy_prof, hOpts, cOptsNoLeg, "", opts.out_dir, "run_avg_vy", "", "runID", "<V_{y}>");
//
//  // Vz - VPD Vz as a function of luminosity
//  TH2D* zdc_dvz = (TH2D*) zdc_vzvpdvz["pre"]->Clone("zdc_dvz");
//  zdc_dvz->Add(zdc_vzvpdvz["low"]);
//  zdc_dvz->Add(zdc_vzvpdvz["mid"]);
//  zdc_dvz->Add(zdc_vzvpdvz["high"]);
//  TH1D* low_lumi_dvz = zdc_dvz->ProjectionY("_lowdvz", 1, 33);
//  low_lumi_dvz->Scale(1.0/low_lumi_dvz->Integral());
//  TH1D* mid_lumi_dvz = zdc_dvz->ProjectionY("_middvz", 34, 66);
//  mid_lumi_dvz->Scale(1.0/mid_lumi_dvz->Integral());
//  TH1D* high_lumi_dvz = zdc_dvz->ProjectionY("_highdvz", 67, 100);
//  high_lumi_dvz->Scale(1.0/high_lumi_dvz->Integral());
//  std::vector<TH1D*> lumi_dvz_bins{low_lumi_dvz, mid_lumi_dvz, high_lumi_dvz};
//  std::vector<string> lumi_dvz_bin_names{"0-33 kHz", "34-66 kHz", "67-100 kHz"};
//  Overlay1D(lumi_dvz_bins, lumi_dvz_bin_names, hOpts, cOpts, opts.out_dir, "lumidvz", "",
//            "V_{z} - VPD V_{z}", "fraction", "ZDC Rate");
//
//
//  // tracks
//  // ------
//
//  // average pt over time
//  TH2D* runid_pt_full = (TH2D*) runid_pt["pre"]->Clone("runid_pt_clone");
//  runid_pt_full->Add(runid_pt["low"]);
//  runid_pt_full->Add(runid_pt["mid"]);
//  runid_pt_full->Add(runid_pt["high"]);
//  TProfile* pt_prof = (TProfile*) runid_pt_full->ProfileX("runid_pt_profile_x");
//  PrettyPrint1D(pt_prof, hOpts, cOptsNoLeg, "", opts.out_dir, "run_avg_pt", "", "runID", "<p_{T}>");
//
//  // pT spectrum with ratios
//  TH2D* zdc_pt_full = (TH2D*) zdc_pt["pre"]->Clone("zdc_pt_clone");
//  zdc_pt_full->Add(zdc_pt["low"]);
//  zdc_pt_full->Add(zdc_pt["mid"]);
//  zdc_pt_full->Add(zdc_pt["high"]);
//  TH1D* low_lumi_pt = zdc_pt_full->ProjectionY("_lowpt", 1, 20);
//  low_lumi_pt->RebinX(4);
//  low_lumi_pt->Scale(1.0/low_lumi_pt->Integral());
//  TH1D* mid_lumi_pt = zdc_pt_full->ProjectionY("_midpt", 21, 66);
//  mid_lumi_pt->RebinX(4);
//  mid_lumi_pt->Scale(1.0/mid_lumi_pt->Integral());
//  TH1D* high_lumi_pt = zdc_pt_full->ProjectionY("_highpt", 67, 100);
//  high_lumi_pt->RebinX(4);
//  high_lumi_pt->Scale(1.0/high_lumi_pt->Integral());
//  TH1D* y7_pt = (TH1D*) zdc_pt["y7"]->ProjectionY("y7_pt");
//  y7_pt->Scale(1.0/y7_pt->Integral());
//  std::vector<TH1D*> tmp_vec{mid_lumi_pt, high_lumi_pt};
//  std::vector<string> tmp_vec_name{"21-66 kHz", "67-100 kHz"};
//  PrintWithRatios(low_lumi_pt, tmp_vec, "0-20 kHz", tmp_vec_name, hOpts, cOpts, opts.out_dir,
//                  "pt_zdc", "", "p_{T}", "fraction", "ZDC Rate");
//
//  // barrel calorimeter
//  // ------------------
//
//  // average tower E_{T}, y14 only
//  TH2D* y14_e = (TH2D*) tow_towe["pre"]->Clone("y14et");
//  y14_e->Add(tow_towe["low"]);
//  y14_e->Add(tow_towe["mid"]);
//  y14_e->Add(tow_towe["high"]);
//  TProfile* y14_avg_e = (TProfile*) y14_e->ProfileX();
//  PrettyPrint1D(y14_avg_e, hOpts, cOptsNoLeg, "", opts.out_dir, "tow_avg_e", "", "towerID", "<E>");
//
//  // average tower energy as a function of run day
//  TH2D* runid_e = (TH2D*) runid_towe["pre"]->Clone("y14runide");
//  runid_e->Add(runid_towe["low"]);
//  runid_e->Add(runid_towe["mid"]);
//  runid_e->Add(runid_towe["high"]);
//  TProfile* y14_runid_avg_e = (TProfile*) runid_e->ProfileX();
//  PrettyPrint1D(y14_runid_avg_e, hOpts, cOptsNoLeg, "", opts.out_dir, "runid_avg_e", "", "runID", "<E>");
//
//  // eT spectra
//  TH2D* zdc_et = (TH2D*) zdc_towet["pre"]->Clone("y14zdcet");
//  zdc_et->Add(zdc_towet["low"]);
//  zdc_et->Add(zdc_towet["mid"]);
//  zdc_et->Add(zdc_towet["high"]);
//  TH1D* low_lumi_et = zdc_et->ProjectionY("_lowet", 1, 33);
//  low_lumi_et->Scale(1.0/low_lumi_et->Integral());
//  TH1D* mid_lumi_et = zdc_et->ProjectionY("_midet", 34, 66);
//  mid_lumi_et->Scale(1.0/mid_lumi_et->Integral());
//  TH1D* high_lumi_et = zdc_et->ProjectionY("_highet", 67, 100);
//  high_lumi_et->Scale(1.0/high_lumi_et->Integral());
//  std::vector<TH1D*> lumi_et_bins{low_lumi_et, mid_lumi_et, high_lumi_et};
//  std::vector<string> lumi_et_bin_names{"0-33 kHz", "34-66 kHz", "67-100 kHz"};
//  Overlay1D(lumi_et_bins, lumi_et_bin_names, hOpts, cOpts, opts.out_dir, "lumiet", "",
//            "E_{T}", "fraction", "ZDC Rate");
  
  return 0;
}
