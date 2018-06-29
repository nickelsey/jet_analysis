// global_primary.cc

#include "jet_analysis/util/common.hh"
#include "jet_analysis/util/string_util.hh"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TF1.h"

#include "TCanvas.h"

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

// create an nprimary/nglobal cut

DEFINE_string(name, "nprim_nglob", "output name");
DEFINE_string(input, "", "input root file");
DEFINE_string(outdir, "tmp", "output directory");
DEFINE_int32(rebinGlobal, 1, "rebin the nglobal axis");
DEFINE_int32(rebinPrimary, 1, "rebin the nprimary axis");
DEFINE_double(nsigma, 2.0, "number of sigma to cut on");

int main(int argc, char* argv[]) {
  
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);
  
  // check to make sure the input file exists
  if (!boost::filesystem::exists(FLAGS_input)) {
    std::cerr << "input file does not exist: " << FLAGS_input << std::endl;;
    return 1;
  }
  // load input
  TFile in(FLAGS_input.c_str(), "READ");
  TH2D* nprimnglob = (TH2D*) in.Get("nglobnprim");
  nprimnglob->RebinX(FLAGS_rebinGlobal);
  nprimnglob->RebinY(FLAGS_rebinPrimary);
  
  // build output directory if it doesn't exist, using boost::filesystem
  if (FLAGS_outdir.empty())
    FLAGS_outdir = "tmp";
  boost::filesystem::path dir(FLAGS_outdir.c_str());
  boost::filesystem::create_directories(dir);
 
  // create output file from the given directory & name
  string outfile_name = FLAGS_outdir + "/" + FLAGS_name + ".root";
  TFile out(outfile_name.c_str(), "RECREATE");
  TH1D* result = new TH1D("nglobnprimcut", "", nprimnglob->GetNbinsX(), 0, nprimnglob->GetXaxis()->GetXmax());
 
  // build our TF1s and TH1s
  int number_of_bins = nprimnglob->GetNbinsX();
  std::vector<TH1D*> projections(number_of_bins);
  std::vector<TF1*> fits(number_of_bins);
 
  for (int bin = 1; bin <= number_of_bins; ++bin) {
 
    projections.push_back(nprimnglob->ProjectionY(MakeString("proj", bin).c_str(), bin, bin));
    
    if (projections.back()->GetEntries() == 0)
      continue;
    
    TH1D* h = projections.back();
    h->Scale(1.0 / projections.back()->Integral());
    h->GetXaxis()->SetRange(3, h->GetNbinsX());
    
    int last_bin = h->FindLastBinAbove(0);
    double last_bin_center = h->GetBinCenter(last_bin);
    int max_bin = h->GetMaximumBin();
    double max_bin_center = h->GetBinCenter(max_bin);
    
    fits.push_back(new TF1(MakeString("fit", bin).c_str(), "gaus(0)", 0, 1500));
    fits.back()->SetParameter(1, 0.0001);
    fits.back()->SetParameter(2, max_bin_center);
    fits.back()->SetParameter(3, (max_bin_center - last_bin_center)/2.0);
    
    LOG(INFO) << "max: " << max_bin_center;
    h->Fit(fits.back(), "EM", "", max_bin_center, 1500);
    
    result->SetBinContent(bin, fits.back()->GetParameter("Mean"));
    result->SetBinError(bin, fits.back()->GetParameter("Sigma"));
    
//    TCanvas c;
//    h->Draw();
//    fits.back()->SetRange(0, 1500);
//    fits.back()->Draw("SAME");
//    c.SaveAs(MakeString("h", bin, ".pdf").c_str());
  }
  
  out.Write();
  out.Close();
  
  gflags::ShutDownCommandLineFlags();
  return 0;
}
