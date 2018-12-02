#include "jet_analysis/util/common.hh"
#include "jet_analysis/util/arg_helper.hh"
#include "jet_analysis/util/string_util.hh"
#include "jet_analysis/util/root_print_routines.hh"
#include "jet_analysis/centrality/centrality_run14.hh"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TStyle.h"
#include "TMath.h"

#include "TCanvas.h"

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

// use glog and gflags
#include "gflags/gflags.h"
#include "glog/stl_logging.h"
#include "glog/logging.h"


DEFINE_string(inputy7, "run7.root", "run 7 data");
DEFINE_string(inputp16, "p16id.root", "P16id data");
DEFINE_string(inputp17, "p17id.root", "P17id data");
DEFINE_string(output, "pt_compare_out", "output directory");

template<class T>
std::vector<TH1D*> SplitOnXAxis(T* h) {
  std::vector<TH1D*> ret;
  
  for (int i = 1; i <= h->GetXaxis()->GetNbins(); ++i) {
    string name = MakeString(h->GetName(), "_", i-1);
    TH1D* tmp = h->ProjectionY(name.c_str(), i, i);
    ret.push_back(tmp);
  }
  
  return ret;
}


int main(int argc, char* argv[]) {
  
  string usage = "generate efficiency curves from embedding QA";
  
  gflags::SetUsageMessage(usage);
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);
  
  // set drawing preferences for histograms and graphs
  gStyle->SetOptStat(false);
  gStyle->SetOptFit(false);
  gStyle->SetLegendBorderSize(0);
  
  // load root file and create output
  if (!boost::filesystem::exists(FLAGS_inputy7)) {
    std::cerr << "input file does not exist: " << FLAGS_inputy7 << std::endl;;
    return 1;
  }
  
  if (!boost::filesystem::exists(FLAGS_inputp16)) {
    std::cerr << "input file does not exist: " << FLAGS_inputp16 << std::endl;;
    return 1;
  }
  if (!boost::filesystem::exists(FLAGS_inputp17)) {
    std::cerr << "input file does not exist: " << FLAGS_inputp17 << std::endl;;
    return 1;
  }
  
  // build output directory if it doesn't exist, using boost::filesystem
  if (FLAGS_output.empty())
  FLAGS_output = "tmp";
  boost::filesystem::path dir(FLAGS_output.c_str());
  boost::filesystem::create_directories(dir);
 
  TFile input_y7(FLAGS_inputy7.c_str(), "READ");
  TFile input_p16(FLAGS_inputp16.c_str(), "READ");
  TFile input_p17(FLAGS_inputp17.c_str(), "READ");
  
  // load input
  TH1D* p16_cent = (TH1D*) input_p16.Get("centrality");
  p16_cent->SetName("p16cent");
  TH2D* p16_pt = (TH2D*) input_p16.Get("pt");
  p16_pt->SetName("ptp16");
  TH2D* p16_pt_corr = (TH2D*) input_p16.Get("ptcorr");
  p16_pt_corr->SetName("ptcorrp16");
  std::vector<TH1D*> p16_pt_cent_bins = SplitOnXAxis(p16_pt);
  std::vector<TH1D*> p16_pt_corr_cent_bins = SplitOnXAxis(p16_pt_corr);
  for (int i = 1; i <= p16_cent->GetXaxis()->GetNbins(); ++i) {
    p16_pt_cent_bins[i-1]->Scale(1.0 / p16_cent->GetBinContent(i));
    p16_pt_corr_cent_bins[i-1]->Scale(1.0 / p16_cent->GetBinContent(i));
  }
  TProfile2D* p16_effic = (TProfile2D*) input_p16.Get("aveeffic");
  p16_effic->SetName("efficp16");
  std::vector<TH1D*> p16_effic_cent_bins = SplitOnXAxis(p16_effic);
  std::reverse(p16_pt_cent_bins.begin(), p16_pt_cent_bins.end());
  std::reverse(p16_pt_corr_cent_bins.begin(), p16_pt_corr_cent_bins.end());
  std::reverse(p16_effic_cent_bins.begin(), p16_effic_cent_bins.end());
  
  TH1D* p17_cent = (TH1D*) input_p17.Get("centrality");
  p17_cent->SetName("p17cent");
  TH2D* p17_pt = (TH2D*) input_p17.Get("pt");
  p17_pt->SetName("p17pt");
  TH2D* p17_pt_corr = (TH2D*) input_p17.Get("ptcorr");
  p17_pt_corr->SetName("p17ptcorr");
  std::vector<TH1D*> p17_pt_cent_bins = SplitOnXAxis(p17_pt);
  std::vector<TH1D*> p17_pt_corr_cent_bins = SplitOnXAxis(p17_pt_corr);
  for (int i = 1; i <= p17_cent->GetXaxis()->GetNbins(); ++i) {
    p17_pt_cent_bins[i-1]->Scale(1.0 / p17_cent->GetBinContent(i));
    p17_pt_corr_cent_bins[i-1]->Scale(1.0 / p17_cent->GetBinContent(i));
  }
  TProfile2D* p17_effic = (TProfile2D*) input_p17.Get("aveeffic");
  p17_effic->SetName("p17effic");
  std::vector<TH1D*> p17_effic_cent_bins = SplitOnXAxis(p17_effic);
  
  TH1D* run7_cent = (TH1D*) input_y7.Get("centrality");
  run7_cent->SetName("run7cent");
  TH2D* run7_pt = (TH2D*) input_y7.Get("pt");
  run7_pt->SetName("run7pt");
  TH2D* run7_pt_corr = (TH2D*) input_y7.Get("ptcorr");
  run7_pt_corr->SetName("run7ptcorr");
  std::vector<TH1D*> run7_pt_cent_bins = SplitOnXAxis(run7_pt);
  std::vector<TH1D*> run7_pt_corr_cent_bins = SplitOnXAxis(run7_pt_corr);
//  run7_pt_corr_cent_bins[0]->Scale(1.0/ run7_cent->GetBinContent(1));
//  run7_pt_corr_cent_bins[1]->Scale(1.0/ run7_cent->GetBinContent(2));
//  run7_pt_corr_cent_bins[2]->Scale(1.0/ run7_cent->GetBinContent(3));
  for (int i = 1; i <= run7_cent->GetXaxis()->GetNbins(); ++i) {
    run7_pt_cent_bins[i-1]->Scale(1.0 / run7_cent->GetBinContent(i));
    run7_pt_corr_cent_bins[i-1]->Scale(1.0 / run7_cent->GetBinContent(i));
  }
  std::vector<TProfile*> run7_eff_cent_bins;
  for (int i = 0; i < 9; ++i) {
    LOG(INFO) << input_y7.Get(MakeString("eff", i).c_str());
    run7_eff_cent_bins.push_back((TProfile*) input_y7.Get(MakeString("eff", i).c_str()));
  }
  
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
  
  
  // do efficiencies
  
  std::vector<TH1D*> print_eff{run7_eff_cent_bins[0], p16_effic_cent_bins[0], p17_effic_cent_bins[0]};
  std::vector<string> name_list{"Run 7", "P16id", "P17id"};
  
  Overlay1D(print_eff, name_list, hopts, cOptsBottomLeg, FLAGS_output, "efficiencies", "", "p_{T}", "efficiency", "0-5% central");
  
  // overlay pTs
  PrintWithRatio(p16_pt_corr_cent_bins[0], p17_pt_corr_cent_bins[0], "P16id", "P17id", hopts, coptslogy, FLAGS_output, "p16_p17",
                 "", "p_{T}", "1/N_{events}dN/dp_{T}", "0-5% central");
  
  PrintWithRatio(run7_pt_corr_cent_bins[0], p17_pt_corr_cent_bins[0], "Run 7", "P17id", hopts, coptslogy, FLAGS_output, "run7_p17",
                 "", "p_{T}", "1/N_{events}dN/dp_{T}", "0-5% central");
  
  PrintWithRatio(run7_pt_corr_cent_bins[0], p16_pt_corr_cent_bins[0], "Run 7", "P16id", hopts, coptslogy, FLAGS_output, "run7_p16",
                 "", "p_{T}", "1/N_{events}dN/dp_{T}", "0-5% central");
  
  LOG(INFO) << "run7: " << run7_pt_corr_cent_bins[0]->Integral();
  LOG(INFO) << "run 14 P17id: " << p17_pt_corr_cent_bins[0]->Integral();
  LOG(INFO) << "before correction";
  LOG(INFO) << "run7: " << run7_pt_cent_bins[0]->Integral();
  LOG(INFO) << "run 14 P17id: " << p17_pt_cent_bins[0]->Integral();
  return 0;
}

//emplate<typename H>
//void Overlay1D(const std::vector<H*>& h,
//               std::vector<std::string> hist_titles,
//               histogramOpts hopts,
//               canvasOpts copts,
//               std::string output_loc,
//               std::string output_name,
//               std::string canvas_title,
//               std::string x_axis_label,
//               std::string y_axis_label,
//               std::string legend_title = "",
//               bool find_good_range = false)

// print histograms & their ratios
//template <typename H>
//void PrintWithRatio(H* h1,
//                    H* h2,
//                    std::string h1_title,
//                    std::string h2_title,
//                    histogramOpts hopts,
//                    canvasOpts copts,
//                    std::string output_loc,
//                    std::string output_name,
//                    std::string canvas_title,
//                    std::string x_axis_label,
//                    std::string y_axis_label,
//                    std::string legend_title = "")
