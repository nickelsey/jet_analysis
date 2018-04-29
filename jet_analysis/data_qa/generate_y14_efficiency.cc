// generate_y14_efficiency.cc
// generate efficiency curves from embedding
// data: two sets of histograms are given,
// mc tracks and reconstructed tracks in bins
// of luminosity, centrality, pt, eta and phi
// this division of reco / mc gives effective
// efficiency

#include "jet_analysis/util/arg_helper.hh"
#include "jet_analysis/util/string_util.hh"
#include "jet_analysis/util/root_print_routines.hh"

#include <string>
#include <vector>
#include <iostream>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TStyle.h"

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

using std::string;

struct Options {
  string name         = "y14_effic";/* name for output root file */
  string input        = "";         /* root file, should contain both y6 & y12 data (with different histogram
                                       prefixes) */
  string out_dir      = "tmp";      /* directory to save output in */
  int lumi_bins       = 3;          /* number of bins in luminosity */
  int cent_bins       = 16;         /* number of bins in centrality */
  double max_pt       = 5.0;        /* maximum pt cutoff */
};

int main(int argc, char* argv[]) {
  
  std::vector<string> refcent_string{"0-5%", "5-10%", "10-15%", "15-20%",
    "20-25%", "25-30%", "30-35%", "35-40%", "40-45%", "45-50%", "50-55%",
    "55-60%", "60-65%", "65-70%", "70-75%", "75-80%"};
  std::vector<string> lumi_string{"0-33 khz", "33-66 khz", "66-100 khz"};
  
  // set drawing preferences for histograms and graphs
  gStyle->SetOptStat(false);
  gStyle->SetOptFit(false);
  gStyle->SetLegendBorderSize(0);
 
  // parse command line options
  Options opts;
  for (int i = 1; i < argc; ++i) {
    if (ParseStrFlag(string(argv[i]), "--name", &opts.name) ||
        ParseStrFlag(string(argv[i]), "--input", &opts.input) ||
        ParseStrFlag(string(argv[i]), "--outDir", &opts.out_dir) ||
        ParseIntFlag(string(argv[i]), "--lumiBins", &opts.lumi_bins) ||
        ParseIntFlag(string(argv[i]), "--centBins", &opts.cent_bins) ||
        ParseFloatFlag(string(argv[i]), "--maxPt", &opts.max_pt)) continue;
    std::cerr << "Unknown command line option: " << argv[i] << std::endl;
    return 1;
  }
  
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
  
  // create output file from the given directory, name & id
  string outfile_name = opts.out_dir + "/" + opts.name + ".root";
  TFile out(outfile_name.c_str(), "RECREATE");
  
  // load input file and read in histograms
  TFile in(opts.input.c_str(), "READ");
  std::vector<std::vector<TH3D*>> in_mc;
  std::vector<std::vector<TH3D*>> in_match;
  
  for (int lumi = 0; lumi < opts.lumi_bins; ++lumi) {
    in_mc.push_back(std::vector<TH3D*>());
    in_match.push_back(std::vector<TH3D*>());
    for (int cent = 0; cent < opts.cent_bins; ++cent) {
      
      string mc_name = "mc_lumi_" + std::to_string(lumi) + "_cent_" + std::to_string(cent);
      string match_name = "match_lumi_" + std::to_string(lumi) + "_cent_" + std::to_string(cent);
      
      in_mc[lumi].push_back((TH3D*) in.Get(mc_name.c_str()));
      in_match[lumi].push_back((TH3D*) in.Get(match_name.c_str()));
      
      in_mc[lumi][cent]->GetXaxis()->SetRange(1, 50);
      in_match[lumi][cent]->GetXaxis()->SetRange(1, 50);
    }
  }
  
  
  // get 2D histograms - integrate out phi dependence for now
  std::vector<std::vector<TH2D*>> mc_2d;
  std::vector<std::vector<TH2D*>> match_2d;
  std::vector<std::vector<TH2D*>> eff_curves;
  std::vector<std::vector<TH2D*>> eff_curves_ratio;
  std::vector<std::vector<TH1D*>> eff_curves_1d;
  std::vector<std::vector<TH1D*>> eff_curves_1d_ratio;
  TCanvas c;
  for (int i = 0; i < opts.lumi_bins; ++i) {
    
    mc_2d.push_back(std::vector<TH2D*>());
    match_2d.push_back(std::vector<TH2D*>());
    eff_curves.push_back(std::vector<TH2D*>());
    eff_curves_ratio.push_back(std::vector<TH2D*>());
    eff_curves_1d.push_back(std::vector<TH1D*>());
    eff_curves_1d_ratio.push_back(std::vector<TH1D*>());
    
    for (int j = 0; j < opts.cent_bins; ++j) {
      mc_2d[i].push_back((TH2D*) in_mc[i][j]->Project3D("YX"));
      match_2d[i].push_back((TH2D*) in_match[i][j]->Project3D("YX"));
      
      string eff_name = "efficiency_lumi_" + std::to_string(i) + "_cent_" + std::to_string(j);
      
      eff_curves[i].push_back((TH2D*) match_2d[i][j]->Clone(eff_name.c_str()));
      eff_curves[i][j]->Divide(mc_2d[i][j]);
      eff_curves_ratio[i].push_back((TH2D*) eff_curves[i][j]->Clone());
      
      string title = lumi_string[i] + " " + refcent_string[j] + " central";
      
      eff_curves[i][j]->SetTitle(title.c_str());
      eff_curves[i][j]->GetZaxis()->SetTitle("efficiency");
      eff_curves[i][j]->GetXaxis()->SetRangeUser(0.0, opts.max_pt);
      eff_curves[i][j]->GetZaxis()->SetRangeUser(0.0, 1.05);
      eff_curves[i][j]->Draw("surf1");
      string full_eff_name = opts.out_dir + "/" + eff_name + ".pdf";
      c.SaveAs(full_eff_name.c_str());
      
      eff_curves_1d[i].push_back(match_2d[i][j]->ProjectionX());
      eff_curves_1d[i][j]->Divide(mc_2d[i][j]->ProjectionX());
      eff_curves_1d_ratio[i].push_back((TH1D*) eff_curves_1d[i][j]->Clone());
      if (i > 0) {
        eff_curves_ratio[i][j]->Divide(eff_curves_ratio[0][j]);
        eff_curves_1d_ratio[i][j]->Divide(eff_curves_1d_ratio[0][j]);
        
        full_eff_name = opts.out_dir + "/" + eff_name + "_ratio.pdf";
        eff_curves_ratio[i][j]->Draw("COLZ");
        c.SaveAs(full_eff_name.c_str());
      }
      
    }
  }
  
  histogramOpts hOpts;
  canvasOpts cOpts;
  cOpts.leg_upper_bound = 0.55;
  cOpts.leg_left_bound = 0.9;
  cOpts.leg_lower_bound = 0.20;
  cOpts.leg_right_bound = 0.5;
  
  for (int i = 0; i < eff_curves_1d.size(); ++i) {
    
    std::string canvas_name;
    std::string file_name;
    if (i == 0) {
      canvas_name = "low luminosity";
      file_name = "low_lumi_eff";
    }
    else if (i == 1) {
      canvas_name = "mid luminosity";
      file_name = "mid_lumi_eff";
    }
    else if (i == 2) {
      canvas_name = "high luminosity";
      file_name = "high_lumi_eff";
    }
    else {
      canvas_name = "Centrality";
      file_name = "centrality_eff";
    }
    eff_curves_1d[i][0]->GetYaxis()->SetRangeUser(0.0, 1.05);
    Overlay1D(eff_curves_1d[i], refcent_string, hOpts, cOpts, opts.out_dir, file_name,
              "", "p_{T}", "efficiency", canvas_name);
    
    if (i != 0) {
      if (i == 1) {
        canvas_name = "mid / low";
        file_name = "mid_low_ratio";
      }
      else if (i == 2) {
        canvas_name = "high / low";
        file_name = "high_low_ratio.pdf";
      }
      Overlay1D(eff_curves_1d_ratio[i], refcent_string, hOpts, cOpts, opts.out_dir, file_name,
                "", "p_{T}", "efficiency", canvas_name);
    }
  }
  
  // get errors
  std::vector<std::vector<TH2D*>> errors;
  for ( int i = 0; i < opts.lumi_bins; ++i) {
    errors.push_back(std::vector<TH2D*>());
    for (int j = 0; j < opts.cent_bins; ++j) {
      string err_name = "error_lumi_" + std::to_string(i) + "_cent_" + std::to_string(j);
      errors[i].push_back((TH2D*) eff_curves[i][j]->Clone(err_name.c_str()));
      
      for (int k = 0; k < (errors[i][j]->GetNbinsX() +2) * (errors[i][j]->GetNbinsY()+1); ++k) {
        errors[i][j]->SetBinContent(k, errors[i][j]->GetBinError(k));
      }
      
      string title = lumi_string[i] + " " + refcent_string[j] + " central: error";
      
      err_name = opts.out_dir + "/" + err_name + ".pdf";
      errors[i][j]->SetTitle(title.c_str());
      errors[i][j]->GetZaxis()->SetRangeUser(0.0, 0.1);
      errors[i][j]->Draw("colz");
      c.SaveAs(err_name.c_str());
    }
  }
  
  
  
  // other histograms for QA
  TH3D* ptmatched = (TH3D*) in.Get("mcptvsmatchptvseta");
  ptmatched->GetXaxis()->SetRangeUser(0.0, 5.0);
  ptmatched->GetYaxis()->SetRangeUser(0.0, 6.0);
  ((TH2D*)ptmatched->Project3D("YX"))->Draw("colz");
  c.SetLogz();
  string matched_name = opts.out_dir + "/" + "match_pt.pdf";
  c.SaveAs(matched_name.c_str());
  
  TH2D* mcvsmatch = (TH2D*) in.Get("mcvsmatched");
  mcvsmatch->Draw("colz");
  string mcvsmatch_name = opts.out_dir + "/" + "mcvsmatch.pdf";
  c.SaveAs(mcvsmatch_name.c_str());
  
  TH2D* refzdc = (TH2D*) in.Get("refzdc");
  string refzdc_name = opts.out_dir + "/" + "refzdc.pdf";
  refzdc->Draw("colz");
  c.SaveAs(refzdc_name.c_str());
  c.SetLogz(false);
  
  TH1D* fit = (TH1D*) in.Get("fitpoints");
  string fit_name = opts.out_dir + "/" + "fitpoints.pdf";
  fit->Draw();
  c.SaveAs(fit_name.c_str());
  
  TH1D* dca = (TH1D*) in.Get("dca");
  string dca_name = opts.out_dir + "/" + "dca.pdf";
  dca->GetXaxis()->SetRangeUser(0.0, 2.0);
  dca->Draw();
  c.SetLogy();
  c.SaveAs(dca_name.c_str());
  c.SetLogy(false);
  // write to file
  out.cd();
  for (auto vec : eff_curves) {
    for (auto hist : vec) {
      hist->Write();
    }
  }
  
  out.Close();
  return 0;
}
