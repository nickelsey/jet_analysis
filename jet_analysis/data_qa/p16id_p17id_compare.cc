// p16id_p17id_compare.cc

#include "jet_analysis/util/common.hh"
#include "jet_analysis/util/string_util.hh"
#include "jet_analysis/util/root_print_routines.hh"

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

#include "TCanvas.h"

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

// compare P16id and P17id on an event-by-event level
// from trees produced by the muDsts

DEFINE_string(name, "p16id_p17id", "output name");
DEFINE_string(p16id, "p16id", "input root file for P16id");
DEFINE_string(p17id, "p17id", "input root file for P17id");
DEFINE_string(tree, "tree", "input tree name");
DEFINE_string(outdir, "tmp", "output directory");

struct Event {
  int runid;
  int eventid;
  int nvertices;
  double vpdvz;
  double zdcrate;
  int nglobal;
  int pxl;
  int ist;
  int ssd;
  std::vector<double> vz;
  std::vector<int> nprim;
  std::vector<int> refmult;
  std::vector<double> rank;
};

typedef unordered_map<pair<int, int>, Event, PairHash> event_map;

bool ReadInTree(string filename, string treename, event_map& map) {
  
  TFile input_file(filename.c_str(), "READ");
  if (input_file.IsOpen()) {
    TTreeReader reader(treename.c_str(), &input_file);
    TTreeReaderValue<int> eventid(reader, "eventid");
    TTreeReaderValue<int> runid(reader, "runid");
    TTreeReaderValue<int> nvertices(reader, "nvertices");
    TTreeReaderValue<int> nglobal(reader, "nglobal");
    TTreeReaderValue<int> pxl(reader, "pxl");
    TTreeReaderValue<int> ist(reader, "ist");
    TTreeReaderValue<int> ssd(reader, "ssd");
    TTreeReaderValue<double> vpdvz(reader, "vpdvz");
    TTreeReaderValue<double> zdcrate(reader, "zdcrate");
    TTreeReaderValue<std::vector<double>> vz(reader, "vz");
    TTreeReaderValue<std::vector<double>> rank(reader, "rank");
    TTreeReaderValue<std::vector<int>> nprim(reader, "nprim");
    TTreeReaderValue<std::vector<int>> refmult(reader, "refmult");
    while (reader.Next()) {
      Event evt;
      evt.runid = *runid;
      evt.eventid = *eventid;
      evt.nvertices = *nvertices;
      evt.vpdvz = *vpdvz;
      evt.zdcrate = *zdcrate;
      evt.nglobal = *nglobal;
      evt.pxl = *pxl;
      evt.ist = *ist;
      evt.ssd = *ssd;
      evt.vz = *vz;
      evt.rank = *rank;
      evt.nprim = *nprim;
      evt.refmult = *refmult;
      map[{*runid, *eventid}] = evt;
    }
    
    return true;
  }
  else
    return false;
}

int main(int argc, char* argv[]) {
  
  string usage = "compares P16id and P17id productions on an event-by-event basis";
  
  gflags::SetUsageMessage(usage);
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);
  
  // set drawing preferences for histograms and graphs
  gStyle->SetOptStat(false);
  gStyle->SetOptFit(false);
  gStyle->SetLegendBorderSize(0);
  
  // check to make sure the input file exists
  if (!boost::filesystem::exists(FLAGS_p16id)) {
    std::cerr << "input file does not exist: " << FLAGS_p16id << std::endl;;
    return 1;
  }
  if (!boost::filesystem::exists(FLAGS_p17id)) {
    std::cerr << "input file does not exist: " << FLAGS_p17id << std::endl;;
    return 1;
  }
  
  // build output directory if it doesn't exist, using boost::filesystem
  if (FLAGS_outdir.empty())
    FLAGS_outdir = "tmp";
  boost::filesystem::path dir(FLAGS_outdir.c_str());
  boost::filesystem::create_directories(dir);
  
  event_map p16id;
  event_map p17id;
  LOG(INFO) << "loading P16id tree";
  ReadInTree(FLAGS_p16id, FLAGS_tree, p16id);
  LOG(INFO) << "loading P17id tree";
  ReadInTree(FLAGS_p17id, FLAGS_tree, p17id);
  LOG(INFO) << "successfully loaded in trees";
  
  TH1D* Vz_p16id = new TH1D("p16id_vz", "", 100, -50, 50);
  TH1D* dVz_p16id = new TH1D("p16id_dvz", "", 100, 0, 50);
  TH1D* Vz_p17id = new TH1D("p17id_vz", "", 100, -50, 50);
  TH1D* dVz_p17id = new TH1D("p17id_dvz", "", 100, 0, 50);
  TH1D* dVz = new TH1D("dvz", "", 100, 0, 20);
  TH1D* dVzmatched = new TH1D("dvzmatched", "", 100, 0, 20);
  TH1D* dnprim = new TH1D("dnprim", "", 200, -500, 500);
  TH1D* dnprimmatched = new TH1D("dnprimmatched", "", 200, -500, 500);
  TH1D* dnglobal = new TH1D("dnglobal", "", 200, -500, 500);
  
  
  for (auto& entry : p16id) {
    auto key = entry.first;
    auto p16id_event = entry.second;
    if (p16id_event.pxl != 0 || p16id_event.ist != 0)
      continue;
    if (p17id.find(key) == p17id.end())
      continue;
    auto p17id_event = p17id[key];
    
    Vz_p16id->Fill(p16id_event.vz[0]);
    Vz_p17id->Fill(p17id_event.vz[0]);
    
    dVz_p16id->Fill(fabs(p16id_event.vz[0] - p16id_event.vpdvz));
    dVz_p17id->Fill(fabs(p17id_event.vz[0] - p17id_event.vpdvz));
    
    dVz->Fill(fabs(p16id_event.vz[0] - p17id_event.vz[0]));
    dnprim->Fill(p16id_event.nprim[0] - p17id_event.vz[0]);
    dnglobal->Fill(p16id_event.nglobal - p17id_event.nglobal);
    
    int nprim_p16 = -1;
    double vz_p16 = 0.0;
    int nprim_p17 = -1;
    double vz_p17 = 0.0;
    for (int i = 0; i < p16id_event.vz.size(); ++i) {
      if (fabs(p16id_event.vz[i] - p16id_event.vpdvz) < 3.0) {
        vz_p16 = p16id_event.vz[i];
        nprim_p16 = p16id_event.nprim[i];
        break;
      }
    }
    for (int i = 0; i < p17id_event.vz.size(); ++i) {
      if (fabs(p17id_event.vz[i] - p17id_event.vpdvz) < 3.0) {
        vz_p17 = p17id_event.vz[i];
        nprim_p17 = p17id_event.nprim[i];
        break;
      }
    }
    if (nprim_p16 < 0 || nprim_p17 < 0)
      continue;
    
    dnprimmatched->Fill(nprim_p16 - nprim_p17);
    dVzmatched->Fill(vz_p16 - vz_p17);
    
  }
  
  
  // print results
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
  
  Vz_p16id->Scale(1.0 / Vz_p16id->Integral());
  Vz_p17id->Scale(1.0 / Vz_p17id->Integral());
  dVz_p16id->Scale(1.0 / dVz_p16id->Integral());
  dVz_p17id->Scale(1.0 / dVz_p17id->Integral());
  dVz->Scale(1.0 / dVz->Integral());
  dnprim->Scale(1.0 / dnprim->Integral());
  dnglobal->Scale(1.0 / dnglobal->Integral());
  dnprimmatched->Scale(1.0 / dnprimmatched->Integral());
  dVzmatched->Scale(1.0 / dVzmatched->Integral());
  
  Overlay1D(Vz_p16id, Vz_p17id, "P16id", "P17id", hopts, copts, FLAGS_outdir, "compare_vz",
            "", "V_{z}", "fraction");
  Overlay1D(dVz_p16id, dVz_p17id, "P16id", "P17id", hopts, coptslogy, FLAGS_outdir, "compare_dvz",
            "", "V_{z}", "fraction");
  PrettyPrint1D(dVz, hopts, coptslogy, "dV_{z}", FLAGS_outdir, "dvz", "", "dV_{z}", "fraction");
  PrettyPrint1D(dnprim, hopts, copts, "dN_{primary}", FLAGS_outdir, "dnprim", "", "dN_{primary}", "fraction");
  PrettyPrint1D(dnprimmatched, hopts, coptslogy, "dN_{primary}", FLAGS_outdir, "dnprimmatched", "", "dN_{primary}", "fraction");
  PrettyPrint1D(dnglobal, hopts, coptslogy, "dN_{global}", FLAGS_outdir, "dnglobal", "", "dN_{global}", "fraction");
  PrettyPrint1D(dVzmatched, hopts, coptslogy, "dV_{z}", FLAGS_outdir, "dvzmatched", "", "dV_{z}", "fraction");
  
  LOG(INFO) << "probability to have 10 or more: " << dnglobal->Integral(dnglobal->GetXaxis()->FindBin(10), 200);
  LOG(INFO) << "probability to have 100 or more: " << dnglobal->Integral(dnglobal->GetXaxis()->FindBin(100), 200);
  
  gflags::ShutDownCommandLineFlags();
  return 0;
}

//void Overlay1D(H* h1,
//               H* h2,
//               std::string h1_title,
//               std::string h2_title,
//               histogramOpts hopts,
//               canvasOpts copts,
//               std::string output_loc,
//               std::string output_name,
//               std::string canvas_title,
//               std::string x_axis_label,
//               std::string y_axis_label,
//               std::string legend_title = "") 

//void PrettyPrint1D(H* h,
//                   histogramOpts hopts,
//                   canvasOpts copts,
//                   std::string hist_title,
//                   std::string output_loc,
//                   std::string output_name,
//                   std::string canvas_title,
//                   std::string x_axis_label,
//                   std::string y_axis_label,
//                   std::string legend_title = "") 
