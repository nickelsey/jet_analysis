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

// compare multiple productions of y14 on an
// event-by-event basis

DEFINE_string(input, "eff_input.root", "root file with embedding QA");
DEFINE_string(name, "efficiency_curves_scaled", "ouput root file name");
DEFINE_string(output, "tmp", "output directory");
DEFINE_double(dcaCut, 1.0, "dca boundary");

using std::string;

std::vector<TH1D*> SplitOnYAxis(TH2D* h) {
  std::vector<TH1D*> ret;
  
  for (int i = 1; i <= h->GetXaxis()->GetNbins(); ++i) {
    string name = MakeString(h->GetName(), "_", i-1);
    TH1D* tmp = h->ProjectionY(name.c_str(), i, i);
    ret.push_back(tmp);
  }
  
  return ret;
}

std::vector<double> ScaleFactor(TH3D* data, TH3D* sim) {
  int centBins = data->GetXaxis()->GetNbins();
  
  std::vector<double> ret(centBins);
  for (int i = 0; i < centBins; ++i) {
    int cent_bin = i+1;
    int pt_bin_low = 1;
    int pt_bin_high = -1;
    string data_name =  MakeString("data_", i);
    string sim_name =  MakeString("sim_", i);
    TH1D* data_tmp = (TH1D*) data->ProjectionZ(data_name.c_str(), cent_bin, cent_bin, pt_bin_low, pt_bin_high);
    TH1D* sim_tmp = (TH1D*) sim->ProjectionZ(sim_name.c_str(), cent_bin, cent_bin, pt_bin_low, pt_bin_high);
    data_tmp->Scale(1.0 / data_tmp->Integral());
    sim_tmp->Scale(1.0 / sim_tmp->Integral());
    double data_int = data_tmp->Integral(data_tmp->FindBin(FLAGS_dcaCut), data_tmp->GetXaxis()->GetNbins());
    double sim_int = sim_tmp->Integral(sim_tmp->FindBin(FLAGS_dcaCut), sim_tmp->GetXaxis()->GetNbins());
    
    ret[i] = sim_int - data_int;
  }
  return ret;
}

std::vector<std::vector<double>> ScaleFactorPtDep(TH3D* data, TH3D* sim) {
  int centBins = data->GetXaxis()->GetNbins();
  int nbins = 20;
  double range_low = 0;
  double range_high = 5.0;
  double dpt = (range_high - range_low) / nbins;
  std::vector<std::vector<double>> ret(centBins, std::vector<double>(nbins));
  for (int i = 0; i < centBins; ++i) {
    data->GetXaxis()->SetRange(i+1, i+1);
    sim->GetXaxis()->SetRange(i+1, i+1);
    for (int j = 0; j < nbins; ++i) {
      
      int cent_bin = i+1;
      double pt_low = range_low + (j * dpt);
      double pt_high = range_low + ((j + 1) * dpt);
      int pt_bin_low = data->GetYaxis()->FindBin(pt_low);
      int pt_bin_high = data->GetYaxis()->FindBin(pt_high);
      
      string data_name =  MakeString("data_", i, "_", j);
      string sim_name =  MakeString("sim_", i, "_", j);
      
      TH1D* data_tmp = (TH1D*) data->ProjectionZ(data_name.c_str(), cent_bin, cent_bin, pt_bin_low, pt_bin_high);
      TH1D* sim_tmp = (TH1D*) sim->ProjectionZ(sim_name.c_str(), cent_bin, cent_bin, pt_bin_low, pt_bin_high);
      data_tmp->Scale(1.0 / data_tmp->Integral());
      sim_tmp->Scale(1.0 / sim_tmp->Integral());
      
      double data_int = data_tmp->Integral(data_tmp->FindBin(FLAGS_dcaCut), data_tmp->GetXaxis()->GetNbins());
      double sim_int = sim_tmp->Integral(sim_tmp->FindBin(FLAGS_dcaCut), sim_tmp->GetXaxis()->GetNbins());
      
      ret[i][j] = sim_int - data_int;
      
    }
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
  if (!boost::filesystem::exists(FLAGS_input)) {
    std::cerr << "input file does not exist: " << FLAGS_input << std::endl;;
    return 1;
  }
  
  // build output directory if it doesn't exist, using boost::filesystem
  if (FLAGS_output.empty())
  FLAGS_output = "tmp";
  boost::filesystem::path dir(FLAGS_output.c_str());
  boost::filesystem::create_directories(dir);
  
  TFile input(FLAGS_input.c_str(), "READ");
  TH2D* mc = (TH2D*) input.Get("mctracks");
  TH2D* rc = (TH2D*) input.Get("recotracks");
  
  rc->Divide(mc);
  
  std::vector<TH1D*> centrality_curves = SplitOnYAxis(rc);
  
  std::string outname = FLAGS_output + "/" + FLAGS_name + ".root";
  TFile out(outname.c_str(), "RECREATE");
  
  for (int i = 0; i < centrality_curves.size(); ++i) {
    std::string name = "cent_bin_" + std::to_string(i);
    centrality_curves[i]->SetName(name.c_str());
    centrality_curves[i]->Write();
  }
  
  out.Close();
  return 0;
}

