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
  bool   use_y7      = true;     /* whether or not to include y7 comparisons */
  bool   use_high    = true;     /* flag to include high luminosity y14 data */
  bool   use_mid     = true;     /* flag to include mid luminosity y14 data */
  bool   use_low     = true;     /* flag to include low luminosity y14 data */
  bool   use_pre     = true;     /* flag to include AuAu_200_production_2014 data */
};

int main(int argc, char* argv[]) {
  
  // set drawing preferences for histograms and graphs
  gStyle->SetOptStat(false);
  gStyle->SetOptFit(false);
  gStyle->SetLegendBorderSize(0);
  
  std::vector<string> prefixes{"y7", "pre", "low", "mid", "high"};
  
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
  std::unordered_map<string, TH2D*> zdc_refmult;
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
  std::unordered_map<string, TH2D*> runid_towe;
  std::unordered_map<string, TH2D*> runid_towet;
  std::unordered_map<string, TH2D*> runid_towadc;
  std::unordered_map<string, TH2D*> zdc_towe;
  std::unordered_map<string, TH2D*> zdc_towet;
  std::unordered_map<string, TH2D*> zdc_towadc;
  std::unordered_map<string, TH2D*> tow_towe;
  std::unordered_map<string, TH2D*> tow_towet;
  std::unordered_map<string, TH2D*> tow_towadc;
  
  // read histograms in from input
  for (auto prefix : prefixes) {
    runid_refmult.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "runidrefmult").c_str())});
    runid_grefmult.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "runidgrefmult").c_str())});
    zdc_refmult.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "zdcrefmult").c_str())});
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
    runid_px.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "runidpx").c_str())});
    runid_py.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "runidpy").c_str())});
    runid_pz.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "runidpz").c_str())});
    runid_pt.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "runidpz").c_str())});
    zdc_px.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "zdcpx").c_str())});
    zdc_py.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "zdcpy").c_str())});
    zdc_pz.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "zdcpz").c_str())});
    zdc_pt.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "zdcpt").c_str())});
    runid_dca.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "runiddca").c_str())});
    runid_fit.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "runidfit").c_str())});
    runid_eta.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "runideta").c_str())});
    runid_phi.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "runidphi").c_str())});
    runid_towe.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "runidtowe").c_str())});
    runid_towet.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "runidtowet").c_str())});
    runid_towadc.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "runidtowadc").c_str())});
    zdc_towe.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "zdctowe").c_str())});
    zdc_towet.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "zdctowet").c_str())});
    zdc_towadc.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "zdctowadc").c_str())});
    tow_towe.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "towtowe").c_str())});
    tow_towet.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "towtowet").c_str())});
    tow_towadc.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "towtowadc").c_str())});
  }
  
  
  // printing routines
  // -----------------
  
  // event level
  // -----------
  TH2D* lumi_refmult = (TH2D*) zdc_refmult["pre"]->Clone("zdc_refmult_dist");
  lumi_refmult->Add(zdc_refmult["low"]);
  lumi_refmult->Add(zdc_refmult["mid"]);
  lumi_refmult->Add(zdc_refmult["high"]);
  TH1D* low_lumi_refmult = (TH1D*) lumi_refmult->ProjectionY("_low", 1, 33);
  low_lumi_refmult->Scale(1.0/low_lumi_refmult->Integral());
  TH1D* mid_lumi_refmult = (TH1D*) lumi_refmult->ProjectionY("_mid", 34, 66);
  mid_lumi_refmult->Scale(1.0/mid_lumi_refmult->Integral());
  TH1D* high_lumi_refmult = (TH1D*) lumi_refmult->ProjectionY("_high", 67, 100);
  high_lumi_refmult->Scale(1.0/high_lumi_refmult->Integral());
  std::vector<TH1D*> lumi_refmult_bins{low_lumi_refmult, mid_lumi_refmult, high_lumi_refmult};
  std::vector<string> lumi_refmult_bin_names{"0-33 kHz", "34-66 kHz", "67-100 kHz"};
  Overlay1D(lumi_refmult_bins, lumi_refmult_bin_names, opts.out_dir, "lumirefmult", "",
            "refMult", "fraction", false, true, true, "ZDC Rate");
  
  // vz, vx, vy
  TH2D* runid_vz_full = (TH2D*) runid_vz["pre"]->Clone("vz_full_clone");
  runid_vz_full->Add(runid_vz["low"]);
  runid_vz_full->Add(runid_vz["mid"]);
  runid_vz_full->Add(runid_vz["high"]);
  TProfile* vz_prof = (TProfile*) runid_vz_full->ProfileX("vz_runid_profile");
  PrettyPrint1D(vz_prof, "", opts.out_dir, "run_avg_vz", "", "runID", "<V_{z}>",
                false, false, false, "");
  
  TH2D* runid_vx_full = (TH2D*) runid_vx["pre"]->Clone("vx_full_clone");
  runid_vx_full->Add(runid_vx["low"]);
  runid_vx_full->Add(runid_vx["mid"]);
  runid_vx_full->Add(runid_vx["high"]);
  TProfile* vx_prof = (TProfile*) runid_vx_full->ProfileX("vx_runid_profile");
  PrettyPrint1D(vx_prof, "", opts.out_dir, "run_avg_vx", "", "runID", "<V_{x}>",
                false, false, false, "");
  
  TH2D* runid_vy_full = (TH2D*) runid_vy["pre"]->Clone("vy_full_clone");
  runid_vy_full->Add(runid_vy["low"]);
  runid_vy_full->Add(runid_vy["mid"]);
  runid_vy_full->Add(runid_vy["high"]);
  TProfile* vy_prof = (TProfile*) runid_vy_full->ProfileX("vy_runid_profile");
  PrettyPrint1D(vy_prof, "", opts.out_dir, "run_avg_vy", "", "runID", "<V_{y}>",
                false, false, false, "");
  

  // barrel calorimeter
  // ------------------
  
  // average tower E_{T}, y14 only
  TH2D* y14_e = (TH2D*) tow_towe["pre"]->Clone("y14et");
  y14_e->Add(tow_towe["low"]);
  y14_e->Add(tow_towe["mid"]);
  y14_e->Add(tow_towe["high"]);
  TProfile* y14_avg_e = (TProfile*) y14_e->ProfileX();
  double y14_e_mean;
  double y14_e_RMS;
  
  std::string y14_avg_e_out = opts.out_dir + "/average_e.pdf";
  MeanStdY(y14_avg_e, y14_e_mean, y14_e_RMS);
  PrintWithBounds(y14_avg_e, y14_e_mean, y14_e_RMS, 3, y14_avg_e_out);
  
  
  // total # of tower hits integrated over run period
  
  TH1D* y14_tower_count = y14_e->ProjectionX();
  double y14_count_mean;
  double y14_count_RMS;
  
  std::string y14_tower_count_out = opts.out_dir + "/tower_count.pdf";
  MeanStdY(y14_tower_count, y14_count_mean, y14_count_RMS);
  PrintWithBounds(y14_tower_count, y14_count_mean, y14_count_RMS, 3, y14_tower_count_out);
  
  
  
  
  return 0;
}
