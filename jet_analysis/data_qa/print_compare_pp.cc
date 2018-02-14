// print_compare_pp.cc

// script to print comparisons between y6 and y12 data

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
  string input       = "";        /* root file, should contain both y6 & y12 data (with different histogram
                                     prefixes) */
  string out_dir      = "tmp";    /* directory to save output in */
  string hist_prefix1 = "";       /* histogram prefix for input1 from comprehensive_data_qa */
  string hist_prefix2 = "";       /* histogram prefix for input2 from comprehensive_data_qa */
  string name1        = "y6";     /* string to identify input 1 in legends and labels */
  string name2        = "y12";    /* string to identify input 2 in legends and labels */
};

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
  
  Options opts;
  for (int i = 1; i < argc; ++i) {
    if (ParseStrFlag(string(argv[i]), "--input", &opts.input) ||
        ParseStrFlag(string(argv[i]), "--outDir", &opts.out_dir) ||
        ParseStrFlag(string(argv[i]), "--histPrefix1", &opts.hist_prefix1) ||
        ParseStrFlag(string(argv[i]), "--histPrefix2", &opts.hist_prefix2) ||
        ParseStrFlag(string(argv[i]), "--name1", &opts.name1) ||
        ParseStrFlag(string(argv[i]), "--name2", &opts.name2)) continue;
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
  
  // make a vector of the histogram prefixes
  std::vector<string> prefixes{opts.hist_prefix1, opts.hist_prefix2};
  
  // we will store histograms in a dictionary for
  // easy lookup
  std::unordered_map<string, TH2D*> runid_refmult;
  std::unordered_map<string, TH2D*> runid_grefmult;
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
    zdc_refmult.insert({prefix, (TH2D*) input.Get(MakeString(prefix, "bbcrefmult").c_str())});
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

  // zdc rates
  TH2D* zdc_refmult_1 = (TH2D*) zdc_refmult[opts.hist_prefix1]->Clone("zdc_refmult_clone_1");
  TH2D* zdc_refmult_2 = (TH2D*) zdc_refmult[opts.hist_prefix2]->Clone("zdc_refmult_clone_2");
  
  TH1D* zdc_rate_1 = (TH1D*) zdc_refmult_1->ProjectionX();
  TH1D* zdc_rate_2 = (TH1D*) zdc_refmult_2->ProjectionX();
  
  zdc_rate_1->Scale(1.0/zdc_rate_1->Integral());
  zdc_rate_1->GetXaxis()->SetRange(1, zdc_rate_1->FindLastBinAbove(0));
  zdc_rate_2->Scale(1.0/zdc_rate_2->Integral());
  
  Overlay1D(zdc_rate_1, zdc_rate_2, opts.name1, opts.name2, hOpts, cOpts,
            opts.out_dir, "zdcrate", "", "ZDC Rate [kHz]", "fraction");
  
  // zdc as a function of runid
  TH2D* runid_zdc_1 = (TH2D*) runid_zdc[opts.hist_prefix1]->Clone("runid_zdc_clone_1");
  TH2D* runid_zdc_2 = (TH2D*) runid_zdc[opts.hist_prefix2]->Clone("runid_zdc_clone_2");
  
  TProfile* runid_zdc_profile_1 = (TProfile*) runid_zdc_1->ProfileX("runid_zdc_clone_1_profx",
                                                                    1, -1, "s");
  TProfile* runid_zdc_profile_2 = (TProfile*) runid_zdc_2->ProfileX("runid_zdc_clone_2_profx",
                                                                    1, -1, "s");
  
  PrettyPrint1D(runid_zdc_profile_1, hOpts, cOptsNoLeg, "", opts.out_dir,
                MakeString("run_avg_zdc_", opts.name1).c_str(), "", "runID", "<ZDC Rate [kHz]>");
  PrettyPrint1D(runid_zdc_profile_2, hOpts, cOptsNoLeg, "", opts.out_dir,
                MakeString("run_avg_zdc_", opts.name2).c_str(), "", "runID", "<ZDC Rate [kHz]>");
  
  // refmult
  TH1D* refmult_1 = (TH1D*) zdc_refmult_1->ProjectionY();
  TH1D* refmult_2 = (TH1D*) zdc_refmult_2->ProjectionY();
  refmult_1->Scale(1.0/refmult_1->Integral());
  refmult_1->GetXaxis()->SetRange(1, refmult_1->FindLastBinAbove(0) + 10);
  refmult_2->Scale(1.0/refmult_2->Integral());
  
  Overlay1D(refmult_1, refmult_2, opts.name1, opts.name2, hOpts, cOptsLogy,
            opts.out_dir, "refmult", "", "refMult", "fraction");
  
  // refmult as a function of runID
  TH2D* runid_refmult_1 = (TH2D*) runid_refmult[opts.hist_prefix1]->Clone("runid_refmult_clone_1");
  TH2D* runid_refmult_2 = (TH2D*) runid_refmult[opts.hist_prefix2]->Clone("runid_refmult_clone_2");
  
  TProfile* runid_refmult_profile_1 = (TProfile*) runid_refmult_1->ProfileX("runid_refmult_clone_1_profx",
                                                                            1, -1, "s");
  TProfile* runid_refmult_profile_2 = (TProfile*) runid_refmult_2->ProfileX("runid_refmult_clone_2_profx",
                                                                            1, -1, "s");
  PrettyPrint1D(runid_refmult_profile_1, hOpts, cOptsNoLeg, "", opts.out_dir,
                MakeString("run_avg_refmult_", opts.name1).c_str(), "", "runID", "<refmult>");
  PrettyPrint1D(runid_refmult_profile_2, hOpts, cOptsNoLeg, "", opts.out_dir,
                MakeString("run_avg_refmult_", opts.name2).c_str(), "", "runID", "<refmult>");
  
  
  // Vz & avg vz
  TH2D* runid_vz_1 = (TH2D*) runid_vz[opts.hist_prefix1]->Clone("runid_vz_clone_1");
  TH2D* runid_vz_2 = (TH2D*) runid_vz[opts.hist_prefix2]->Clone("runid_vz_clone_2");
  
  TProfile* runid_vz_profile_1 = (TProfile*) runid_vz_1->ProfileX("runid_vz_clone_1_profx", 1, -1, "s");
  TProfile* runid_vz_profile_2 = (TProfile*) runid_vz_2->ProfileX("runid_vz_clone_2_profx", 1, -1, "s");
  
  PrettyPrint1D(runid_vz_profile_1, hOpts, cOptsNoLeg, "", opts.out_dir,
                MakeString("runid_avg_vz_", opts.name1).c_str(), "", "runID", "<V_{z}>");
  PrettyPrint1D(runid_vz_profile_2, hOpts, cOptsNoLeg, "", opts.out_dir,
                MakeString("runid_avg_vz_", opts.name2).c_str(), "", "runID", "<V_{z}>");
  
  TH1D* vz_1 = (TH1D*) runid_vz_1->ProjectionY();
  TH1D* vz_2 = (TH1D*) runid_vz_2->ProjectionY();
  
  vz_1->Scale(1.0 / vz_1->Integral());
  vz_2->Scale(1.0 / vz_2->Integral());
  
  Overlay1D(vz_1, vz_2, opts.name1, opts.name2, hOpts, cOptsBottomLeg, opts.out_dir,
            "vz", "", "V_{z}", "fraction");
  
  // Vx and Vy
  TH2D* vx_vy_1 = (TH2D*) vx_vy[opts.hist_prefix1]->Clone("vx_vy_clone_1");
  TH2D* vx_vy_2 = (TH2D*) vx_vy[opts.hist_prefix2]->Clone("vx_vy_clone_2");
  
  TH1D* vx_1 = vx_vy_1->ProjectionX();
  TH1D* vx_2 = vx_vy_2->ProjectionX();
  vx_1->Scale(1.0 / vx_1->Integral());
  vx_2->Scale(1.0 / vx_2->Integral());
  
  TH1D* vy_1 = vx_vy_1->ProjectionY();
  TH1D* vy_2 = vx_vy_2->ProjectionY();
  vy_1->Scale(1.0 / vy_1->Integral());
  vy_2->Scale(1.0 / vy_2->Integral());
  
  Overlay1D(vx_1, vx_2, opts.name1, opts.name2, hOpts, cOpts, opts.out_dir,
            "vx", "", "V_{x}", "fraction");
  Overlay1D(vy_1, vy_2, opts.name1, opts.name2, hOpts, cOpts, opts.out_dir,
            "vy", "", "V_{y}", "fraction");
  
  // pT
  TH2D* runid_pt_1 = (TH2D*) runid_pt[opts.hist_prefix1]->Clone("runid_pt_clone_1");
  TH2D* runid_pt_2 = (TH2D*) runid_pt[opts.hist_prefix2]->Clone("runid_pt_clone_2");
  
  TProfile* runid_pt_profile_1 = (TProfile*) runid_pt_1->ProfileX("runid_pt_clone_1_profx",
                                                                  1, -1, "s");
  TProfile* runid_pt_profile_2 = (TProfile*) runid_pt_2->ProfileX("runid_pt_clone_2_profx",
                                                                  1, -1, "s");
  
  PrettyPrint1D(runid_pt_profile_1, hOpts, cOptsNoLeg, "", opts.out_dir,
                MakeString("pt_", opts.name1).c_str(), "", "runID", "<p_{T}>");
  PrettyPrint1D(runid_pt_profile_2, hOpts, cOptsNoLeg, "", opts.out_dir,
                MakeString("pt_", opts.name2).c_str(), "", "runID", "<p_{T}>");
  
  TH1D* pt_1 = (TH1D*) runid_pt_1->ProjectionY();
  TH1D* pt_2 = (TH1D*) runid_pt_2->ProjectionY();
  pt_1->RebinX(2);
  pt_2->RebinX(2);
  pt_1->Scale(1.0 / pt_1->Integral(2, -1));
  pt_2->Scale(1.0 / pt_2->Integral(2, -1));
  
  Overlay1D(pt_1, pt_2, opts.name1, opts.name2, hOpts, cOptsLogy, opts.out_dir,
            "pt", "", "p_{T}", "fraction");
  PrintWithRatio(pt_1, pt_2, opts.name1, opts.name2, hOpts, cOptsLogy, opts.out_dir,
                 "pt_ratio", "", "p_{T}", "fraction");

  
  // px
  TH2D* runid_px_1 = (TH2D*) runid_px[opts.hist_prefix1]->Clone("runid_px_clone_1");
  TH2D* runid_px_2 = (TH2D*) runid_px[opts.hist_prefix2]->Clone("runid_px_clone_2");
  
  TProfile* runid_px_profile_1 = (TProfile*) runid_px_1->ProfileX("runid_px_clone_1_profx",
                                                                  1, -1, "s");
  TProfile* runid_px_profile_2 = (TProfile*) runid_px_2->ProfileX("runid_px_clone_2_profx",
                                                                  1, -1, "s");
  PrettyPrint1D(runid_px_profile_1, hOpts, cOptsNoLeg, "", opts.out_dir,
                MakeString("px_", opts.name1).c_str(), "", "runID", "<p_{x}>");
  PrettyPrint1D(runid_px_profile_2, hOpts, cOptsNoLeg, "", opts.out_dir,
                MakeString("px_", opts.name2).c_str(), "", "runID", "<p_{x}>");
  
  TH1D* px_1 = (TH1D*) runid_px_1->ProjectionY();
  TH1D* px_2 = (TH1D*) runid_px_2->ProjectionY();
  px_1->RebinX(2);
  px_2->RebinX(2);
  px_1->Scale(1.0 / px_1->Integral(2, -1));
  px_2->Scale(1.0 / px_2->Integral(2, -1));
  
  Overlay1D(px_1, px_2, opts.name1, opts.name2, hOpts, cOptsLogy, opts.out_dir,
            "px", "", "p_{x}", "fraction");
  PrintWithRatio(px_1, px_2, opts.name1, opts.name2, hOpts, cOptsLogy, opts.out_dir,
                 "px_ratio", "", "p_{x}", "fraction");
  
  // py
  TH2D* runid_py_1 = (TH2D*) runid_py[opts.hist_prefix1]->Clone("runid_py_clone_1");
  TH2D* runid_py_2 = (TH2D*) runid_py[opts.hist_prefix2]->Clone("runid_py_clone_2");
  
  TProfile* runid_py_profile_1 = (TProfile*) runid_py_1->ProfileX("runid_py_clone_1_profx",
                                                                  1, -1, "s");
  TProfile* runid_py_profile_2 = (TProfile*) runid_py_2->ProfileX("runid_py_clone_2_profx",
                                                                  1, -1, "s");
  
  PrettyPrint1D(runid_py_profile_1, hOpts, cOptsNoLeg, "", opts.out_dir,
                MakeString("py_", opts.name1).c_str(), "", "runID", "<p_{y}>");
  PrettyPrint1D(runid_py_profile_2, hOpts, cOptsNoLeg, "", opts.out_dir,
                MakeString("py_", opts.name2).c_str(), "", "runID", "<p_{y}>");
  
  TH1D* py_1 = (TH1D*) runid_py_1->ProjectionY();
  TH1D* py_2 = (TH1D*) runid_py_2->ProjectionY();
  py_1->RebinX(2);
  py_2->RebinX(2);
  py_1->Scale(1.0 / py_1->Integral(2, -1));
  py_2->Scale(1.0 / py_2->Integral(2, -1));
  
  Overlay1D(py_1, py_2, opts.name1, opts.name2, hOpts, cOptsLogy, opts.out_dir,
            "py", "", "p_{y}", "fraction");
  PrintWithRatio(py_1, py_2, opts.name1, opts.name2, hOpts, cOptsLogy, opts.out_dir,
                 "py_ratio", "", "p_{y}", "fraction");
  
  // pz
  TH2D* runid_pz_1 = (TH2D*) runid_pz[opts.hist_prefix1]->Clone("runid_pz_clone_1");
  TH2D* runid_pz_2 = (TH2D*) runid_pz[opts.hist_prefix2]->Clone("runid_pz_clone_2");
  
  TProfile* runid_pz_profile_1 = (TProfile*) runid_pz_1->ProfileX("runid_pz_clone_1_profx",
                                                                  1, -1, "s");
  TProfile* runid_pz_profile_2 = (TProfile*) runid_pz_2->ProfileX("runid_pz_clone_2_profx",
                                                                  1, -1, "s");
  
  PrettyPrint1D(runid_pz_profile_1, hOpts, cOptsNoLeg, "", opts.out_dir,
                MakeString("pz_", opts.name1).c_str(), "", "runID", "<p_{z}>");
  PrettyPrint1D(runid_pz_profile_2, hOpts, cOptsNoLeg, "", opts.out_dir,
                MakeString("pz_", opts.name2).c_str(), "", "runID", "<p_{z}>");

  TH1D* pz_1 = (TH1D*) runid_pz_1->ProjectionY();
  TH1D* pz_2 = (TH1D*) runid_pz_2->ProjectionY();
  pz_1->RebinX(2);
  pz_2->RebinX(2);
  pz_1->Scale(1.0 / pz_1->Integral(2, -1));
  pz_2->Scale(1.0 / pz_2->Integral(2, -1));
  
  Overlay1D(pz_1, pz_2, opts.name1, opts.name2, hOpts, cOptsLogy, opts.out_dir,
            "pz", "", "p_{z}", "fraction");
  PrintWithRatio(pz_1, pz_2, opts.name1, opts.name2, hOpts, cOptsLogy, opts.out_dir,
                 "pz_ratio", "", "p_{z}", "fraction");
  
  // et
  TH2D* runid_towet_1 = (TH2D*) runid_towet[opts.hist_prefix1]->Clone("runid_towet_clone_1");
  TH2D* runid_towet_2 = (TH2D*) runid_towet[opts.hist_prefix2]->Clone("runid_towet_clone_2");
  
  TProfile* runid_towet_profile_1 = (TProfile*) runid_towet_1->ProfileX("runid_towet_clone_1_profx",
                                                                        1, -1, "s");
  TProfile* runid_towet_profile_2 = (TProfile*) runid_towet_2->ProfileX("runid_towet_clone_2_profx",
                                                                        1, -1, "s");
  
  PrettyPrint1D(runid_towet_profile_1, hOpts, cOptsNoLeg, "", opts.out_dir,
                MakeString("et_", opts.name1).c_str(), "", "runID", "<E_{T}>");
  PrettyPrint1D(runid_towet_profile_2, hOpts, cOptsNoLeg, "", opts.out_dir,
                MakeString("et_", opts.name2).c_str(), "", "runID", "<E_{T}>");
  
  TH1D* et_1 = (TH1D*) runid_towet_1->ProjectionY();
  TH1D* et_2 = (TH1D*) runid_towet_2->ProjectionY();
  et_1->RebinX(2);
  et_2->RebinX(2);
  et_1->Scale(1.0 / et_1->Integral(2, -1));
  et_2->Scale(1.0 / et_2->Integral(2, -1));
  
  Overlay1D(et_1, et_2, opts.name1, opts.name2, hOpts, cOptsLogy, opts.out_dir,
            "et", "", "E_{T}", "fraction");
  PrintWithRatio(et_1, et_2, opts.name1, opts.name2, hOpts, cOptsLogy, opts.out_dir,
                 "et_ratio", "", "E_{T}", "fraction");
  
  return 0;
};
