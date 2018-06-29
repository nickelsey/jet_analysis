// compare_auau_print.cc

#include "jet_analysis/util/common.hh"
#include "jet_analysis/util/string_util.hh"
#include "jet_analysis/util/root_print_routines.hh"

#include <iostream>
#include <exception>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TStyle.h"
#include "TProfile.h"

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

DEFINE_string(input, "auau_compare", "input directory");
DEFINE_string(outdir, "auau_compare/tmp", "output directory");

enum class DATA {y14, y11, y7, p16};
map<DATA, string> data_types{{DATA::y14, "y14"}, {DATA::y11, "y11"},
  {DATA::y7, "y7"}, {DATA::p16, "p16"}};
enum class CUTS {low, mid, high, full};
map<CUTS, string> cut_types{{CUTS::low, "low"}, {CUTS::mid, "mid"},
  {CUTS::high, "high"}, {CUTS::full, "full"}};

map<pair<DATA, CUTS>, TFile*> LoadAllInputs(string inputdir) {
  map<pair<DATA, CUTS>, TFile*> ret;
  for (auto type : data_types) {
    for (auto cut : cut_types) {
      string file = MakeString(inputdir, "/", type.second, "_", cut.second, ".root");
      // check to make sure the input file exists
      if (!boost::filesystem::exists(file))
        throw std::runtime_error(MakeString("could not load file: ", file).c_str());
      ret[{type.first, cut.first}] = new TFile(file.c_str(), "READ");
    }
  }
  return ret;
}

int main(int argc, char* argv[]) {
  
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);
  
  // set drawing preferences for histograms and graphs
  gStyle->SetOptStat(false);
  gStyle->SetOptFit(false);
  gStyle->SetLegendBorderSize(0);
  
  // build output directory if it doesn't exist, using boost::filesystem
  if (FLAGS_outdir.empty())
    FLAGS_outdir = "tmp";
  boost::filesystem::path dir(FLAGS_outdir.c_str());
  boost::filesystem::create_directories(dir);
  
  auto input_files = LoadAllInputs(FLAGS_input);
  
  map<pair<DATA, CUTS>, TH3D*> tracks;
  map<pair<DATA, CUTS>, TH3D*> tracksvz;
  map<pair<DATA, CUTS>, TH3D*> lumitracks;
  map<pair<DATA, CUTS>, TH3D*> lumitracksvz;
  map<pair<DATA, CUTS>, TH2D*> refmult;
  map<pair<DATA, CUTS>, TH2D*> refmultvz;
  map<pair<DATA, CUTS>, TH2D*> pt;
  map<pair<DATA, CUTS>, TH2D*> nhit;
  map<pair<DATA, CUTS>, TH2D*> nhitpos;
  map<pair<DATA, CUTS>, TH2D*> fitfrac;
  map<pair<DATA, CUTS>, TProfile*> avgnglobal;
  map<pair<DATA, CUTS>, TProfile*> avgnhit;
  map<pair<DATA, CUTS>, TH3D*> nglobaldca;
  
  for (auto type : data_types) {
    for (auto cut : cut_types) {
      tracks[{type.first, cut.first}] = (TH3D*) input_files[{type.first, cut.first}]->Get("tracks");
      tracksvz[{type.first, cut.first}] = (TH3D*) input_files[{type.first, cut.first}]->Get("tracksvz");
      lumitracks[{type.first, cut.first}] = (TH3D*) input_files[{type.first, cut.first}]->Get("lumitracks");
      lumitracksvz[{type.first, cut.first}] = (TH3D*) input_files[{type.first, cut.first}]->Get("lumitracksvz");
      refmult[{type.first, cut.first}] = (TH2D*) input_files[{type.first, cut.first}]->Get("refmult");
      refmultvz[{type.first, cut.first}] = (TH2D*) input_files[{type.first, cut.first}]->Get("refmultvz");
      pt[{type.first, cut.first}] = (TH2D*) input_files[{type.first, cut.first}]->Get("pt");
      nhit[{type.first, cut.first}] = (TH2D*) input_files[{type.first, cut.first}]->Get("nhit");
      nhitpos[{type.first, cut.first}] = (TH2D*) input_files[{type.first, cut.first}]->Get("nhitpos");
      fitfrac[{type.first, cut.first}] = (TH2D*) input_files[{type.first, cut.first}]->Get("fitfrac");
      avgnglobal[{type.first, cut.first}] = (TProfile*) input_files[{type.first, cut.first}]->Get("avgnglobal");
      avgnhit[{type.first, cut.first}] = (TProfile*) input_files[{type.first, cut.first}]->Get("avgnhitvz");
      nglobaldca[{type.first, cut.first}] = (TH3D*) input_files[{type.first, cut.first}]->Get("nglobal_dca");
    }
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

  // first we want to look at the reco refmult and dca distribution for all datasets w/ low cuts
  for (auto data : data_types) {
    Print2DSimple(refmultvz[{data.first, CUTS::low}], hopts, coptslogz, FLAGS_outdir, MakeString("refmultvz_", data.second),
                  "", "refmult", "reco refmult", "colz");
    Print2DSimple(((TH2D*)nglobaldca[{data.first, CUTS::low}]->Project3D("YX")), hopts, coptslogz, FLAGS_outdir, MakeString("nglobaldca_", data.second),
                  "", "nglobal", "DCA [cm]", "colz");
  }
  
  // look at refmult for all 4 datasets
  vector<TH1D*> refmult_data;
  vector<TH1D*> refmult_data_alt;
  vector<string> refmult_data_name;
  for (auto data : data_types) {
    refmult_data.push_back(refmultvz[{data.first, CUTS::low}]->ProjectionX(MakeString("refmult", data.second, "low").c_str()));
    refmult_data_alt.push_back(refmultvz[{data.first, CUTS::low}]->ProjectionY(MakeString("refmultalt", data.second, "low").c_str()));
    refmult_data_name.push_back(data.second);
    refmult_data.back()->Scale(1.0 / refmult_data.back()->Integral());
  }
  Overlay1D(refmult_data, refmult_data_name, hopts, coptslogy, FLAGS_outdir, "refmult",
            "", "refmult", "fraction");
  Overlay1D(refmult_data, refmult_data_name, hopts, coptslogy, FLAGS_outdir, "refmultalt",
            "", "refmult", "fraction");
  
  // we want to see the evolution of nhit as a function of cut quality
  for (auto data : data_types) {
    vector<TH1D*> nhit_tmp;
    vector<TH1D*> fitfrac_tmp;
    vector<TH1D*> nhitpos_tmp;
    vector<string> names_tmp;
    for (auto cut : cut_types) {
      nhit_tmp.push_back(nhit[{data.first, cut.first}]->ProjectionY(MakeString("nhit", data.second, cut.second).c_str(), 1, 1));
      nhit_tmp.back()->Scale(1.0 / nhit_tmp.back()->Integral());
      nhitpos_tmp.push_back(nhitpos[{data.first, cut.first}]->ProjectionY(MakeString("nhitpos", data.second, cut.second).c_str(), 1, 1));
      nhitpos_tmp.back()->Scale(1.0 / nhitpos_tmp.back()->Integral());
      fitfrac_tmp.push_back(fitfrac[{data.first, cut.first}]->ProjectionY(MakeString("fitfrac", data.second, cut.second).c_str(), 1, 1));
      fitfrac_tmp.back()->Scale(1.0 / fitfrac_tmp.back()->Integral());
      names_tmp.push_back(cut.second);
    }
    Overlay1D(nhit_tmp, names_tmp, hopts, cOptsBottomLeftLeg, FLAGS_outdir, MakeString("nhit", data.second), "", "nhit", "fraction");
    Overlay1D(nhitpos_tmp, names_tmp, hopts, cOptsBottomLeftLegLogy, FLAGS_outdir, MakeString("nhitpos", data.second), "", "nhitpos", "fraction");
    Overlay1D(fitfrac_tmp, names_tmp, hopts, cOptsBottomLeftLeg, FLAGS_outdir, MakeString("fitfrac", data.second), "", "fitfrac", "fraction");
  }
  
  // print avg nglobal as a function of nprimary
  for (auto cut : cut_types) {
    vector<TProfile*> avg_nglobal_tmp;
    vector<string> names_tmp;
    for (auto data : data_types) {
      avg_nglobal_tmp.push_back(avgnglobal[{data.first, cut.first}]);
      names_tmp.push_back(data.second);
    }
    Overlay1D(avg_nglobal_tmp, names_tmp, hopts, cOptsBottomLeg, FLAGS_outdir, MakeString("avgnglobal", cut.second), "",
              "nprimary", "<nglobal>");
  }
  
  // print comparison of average nhit
  for (auto cut : cut_types) {
    vector<TProfile*> avg_nhit_tmp;
    vector<string> names_tmp;
    for (auto data : data_types) {
      avg_nhit_tmp.push_back(avgnhit[{data.first, cut.first}]);
      names_tmp.push_back(data.second);
    }
    Overlay1D(avg_nhit_tmp, names_tmp, hopts, copts, FLAGS_outdir, MakeString("avgnhit", cut.second), "",
              "nglobal", "<nhit>");
  }
  
  
  // compare y11 and y14 nhits
  for (auto cut : cut_types) {
    auto h1 = (TH1D*) nhit[{DATA::y11, cut.first}]->ProjectionY(MakeString("tmpnhity11", cut.second).c_str(), 1, 1);
    auto h2 = (TH1D*) nhit[{DATA::y14, cut.first}]->ProjectionY(MakeString("tmpnhity14", cut.second).c_str(), 1, 1);
    h1->Scale(1.0 / h1->Integral());
    h2->Scale(1.0 / h2->Integral());
    //h1->GetXaxis()->SetRange(lowBin, highBin);
    //h2->GetXaxis()->SetRange(lowBin, highBin);
    PrintWithRatio(h1, h2, "y11", "y14", hopts, cOptsBottomLeftLegLogy, FLAGS_outdir, MakeString("nhitratio", cut.second), "", "nhits", "fraction");
  }
  
  // now do average DCA
  vector<pair<int, int>> pt_bins{{1, 10}, {11, 20}, {21, 50}};
  vector<string> pt_bin_names = {"0_1_gev", "1_2_gev", "2_5_gev"};
  for (auto cut : cut_types) {
    for (int i = 0; i < pt_bins.size(); ++i) {
    vector<TProfile*> avg_dca_tmp;
    vector<string> names_tmp;
    for (auto data : data_types) {
      nglobaldca[{data.first, cut.first}]->GetZaxis()->SetRange(pt_bins[i].first, pt_bins[i].second);
      avg_dca_tmp.push_back(((TH2D*)nglobaldca[{data.first, cut.first}]->Project3D("YX"))->ProfileX());
      avg_dca_tmp.back()->SetName(MakeString(avg_dca_tmp.back()->GetName(), cut.second, data.second, pt_bin_names[i]).c_str());
      names_tmp.push_back(data.second);
    }
    if (i == 0)
      avg_dca_tmp[0]->GetYaxis()->SetRangeUser(0.4, 0.8);
    if (i == 1)
      avg_dca_tmp[0]->GetYaxis()->SetRangeUser(0.25, 0.6);
    if (i == 2)
      avg_dca_tmp[0]->GetYaxis()->SetRangeUser(0.2, 0.5);
    Overlay1D(avg_dca_tmp, names_tmp, hopts, cOptsTopLeftLeg, FLAGS_outdir, MakeString("avgdca", cut.second, pt_bin_names[i]), "",
              "nglobal", "<DCA>", pt_bin_names[i], false);
    }
  }
  
  gflags::ShutDownCommandLineFlags();
  return 0;
}

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

//void Overlay1D(const std::vector<H*>& h,
//               std::vector<std::string> hist_titles,
//               histogramOpts hopts,
//               canvasOpts copts,
//               std::string output_loc,
//               std::string output_name,
//               std::string canvas_title,
//               std::string x_axis_label,
//               std::string y_axis_label,
//               std::string legend_title = "")

//void Print2DSimple(H* h,
//                   histogramOpts hopts,
//                   canvasOpts copts,
//                   std::string output_loc,
//                   std::string output_name,
//                   std::string canvas_title,
//                   std::string x_axis_label,
//                   std::string y_axis_label,
//                   std::string opt = "COLZ")
