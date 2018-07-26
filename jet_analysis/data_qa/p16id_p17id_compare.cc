// p16id_p17id_compare.cc

#include "jet_analysis/util/common.hh"
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
  std::vector<double> pt;
  std::vector<double> dca;
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
    TTreeReaderValue<std::vector<double>> dca(reader, "dca");
    TTreeReaderValue<std::vector<double>> pt(reader, "pt");
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
      evt.dca = *dca;
      evt.pt = *pt;
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
 
  TH1D* rank_p16id = new TH1D("p16id_rank", "", 150, -10, 5);
  TH1D* Vz_p16id = new TH1D("p16id_vz", "", 100, -50, 50);
  TH1D* dVz_p16id = new TH1D("p16id_dvz", "", 100, 0, 50);
  TH1D* Vz_p17id = new TH1D("p17id_vz", "", 100, -50, 50);
  TH1D* rank_p17id = new TH1D("p17id_rank", "", 150, -10, 5);
  TH1D* dVz_p17id = new TH1D("p17id_dvz", "", 100, 0, 50);
  TH1D* dVz = new TH1D("dvz", "", 100, 0, 20);
  TH1D* dVzmatched = new TH1D("dvzmatched", "", 100, 0, 20);
  TH1D* dnprim = new TH1D("dnprim", "", 200, -500, 500);
  TH1D* dnprimmatched = new TH1D("dnprimmatched", "", 200, -500, 500);
  TH1D* dnglobal = new TH1D("dnglobal", "", 200, -500, 500);
  
  // some centrality analyses
  TH2D* dVzToNominal_p16 = new TH2D("dvztonom16", "", 8, 0, 8, 100, 0, 50);
  TH2D* nominalDCA_p16   = new TH2D("nominaldca16", "", 8, 0, 8, 100, 0, 5);
  TH2D* nominalPt_p16    = new TH2D("nominalpt16", "", 8, 0, 8, 100, 0, 2);
  TH2D* nominalMult_p16  = new TH2D("nominalmult16", "", 8, 0, 8, 100, 0, 1500);
  TH2D* closestDCA_p16   = new TH2D("closestdca16", "", 8, 0, 8, 100, 0, 5);
  TH2D* closestPt_p16    = new TH2D("closestpt16", "", 8, 0, 8, 100, 0, 2);
  TH2D* closestMult_p16  = new TH2D("closestmult16", "", 8, 0, 8, 100, 0, 1500);
  
  TH2D* dVzToNominal_p17 = new TH2D("dvztonom17", "", 8, 0, 8, 100, 0, 50);
  TH2D* nominalDCA_p17   = new TH2D("nominaldca17", "", 8, 0, 8, 100, 0, 5);
  TH2D* nominalPt_p17    = new TH2D("nominalpt17", "", 8, 0, 8, 100, 0, 2);
  TH2D* nominalMult_p17  = new TH2D("nominalmult17", "", 8, 0, 8, 100, 0, 1500);
  TH2D* closestDCA_p17   = new TH2D("closestdca17", "", 8, 0, 8, 100, 0, 5);
  TH2D* closestPt_p17    = new TH2D("closestpt17", "", 8, 0, 8, 100, 0, 2);
  TH2D* closestMult_p17  = new TH2D("closestmult17", "", 8, 0, 8, 100, 0, 1500);
  
  // create centrality definition
  CentralityRun14 centrality;
  
  for (auto& entry : p16id) {
    auto key = entry.first;
    auto p16id_event = entry.second;
    if (p16id_event.pxl != 0 || p16id_event.ist != 0)
      continue;
    if (p17id.find(key) == p17id.end())
      continue;
    auto p17id_event = p17id[key];

    int id_p16 = -1;
    int nprim_p16 = -1;
    double vz_p16 = 0.0;
    double dca_p16 = 0;
    double pt_p16 = 0;
    int id_p17 = -1;
    int nprim_p17 = -1;
    double vz_p17 = 0.0;
    double dca_p17 = 0;
    double pt_p17 = 0;
    
    for (int i = 0; i < p16id_event.vz.size(); ++i) {
      if (fabs(p16id_event.vz[i] - p16id_event.vpdvz) < 3.0) {
        id_p16 = i;
        vz_p16 = p16id_event.vz[i];
        nprim_p16 = p16id_event.nprim[i];
        dca_p16 = p16id_event.dca[i];
        pt_p16 = p16id_event.pt[i];
        break;
      }
    }
    for (int i = 0; i < p17id_event.vz.size(); ++i) {
      if (fabs(p17id_event.vz[i] - p17id_event.vpdvz) < 3.0) {
        id_p17 = i;
        vz_p17 = p17id_event.vz[i];
        nprim_p17 = p17id_event.nprim[i];
        dca_p17 = p17id_event.dca[i];
        pt_p17 = p17id_event.pt[i];
        break;
      }
    }
    
    if (nprim_p16 < 0 || nprim_p17 < 0)
      continue;

    rank_p16id->Fill(p16id_event.rank[id_p16]);
    rank_p17id->Fill(p17id_event.rank[id_p17]);
   
    Vz_p16id->Fill(p16id_event.vz[id_p16]);
    Vz_p17id->Fill(p17id_event.vz[id_p17]);
    
    dVz_p16id->Fill(fabs(p16id_event.vz[id_p16] - p16id_event.vpdvz));
    dVz_p17id->Fill(fabs(p17id_event.vz[id_p17] - p17id_event.vpdvz));
    
    dVz->Fill(fabs(p16id_event.vz[id_p16] - p17id_event.vz[id_p17]));
    dnprim->Fill(p16id_event.nprim[id_p16] - p17id_event.vz[id_p17]);
    dnglobal->Fill(p16id_event.nglobal - p17id_event.nglobal);
    


    centrality.setEvent(p16id_event.runid, p16id_event.refmult[id_p16], p16id_event.zdcrate, vz_p16);
    int cent16 = centrality.centrality9();
    centrality.setEvent(p17id_event.runid, p17id_event.refmult[id_p17], p17id_event.zdcrate, vz_p17);
    int cent17 = centrality.centrality9();
    
    nominalDCA_p16->Fill(cent16, dca_p16);
    nominalPt_p16->Fill(cent16, pt_p16);
    nominalMult_p16->Fill(cent16, nprim_p16);
    
    nominalDCA_p17->Fill(cent17, dca_p17);
    nominalPt_p17->Fill(cent17, pt_p17);
    nominalMult_p17->Fill(cent17, nprim_p17);
    
    dnprimmatched->Fill(nprim_p16 - nprim_p17);
    dVzmatched->Fill(vz_p16 - vz_p17);
    
    // find the secondary vertices
    int nprim_secondary_p16 = -1;
    int nprim_secondary_p17 = -1;
    double vz_secondary_p16 = 0;
    double vz_secondary_p17 = 0;
    double dca_secondary_p16 = 0;
    double dca_secondary_p17 = 0;
    double pt_secondary_p16 = 0;
    double pt_secondary_p17 = 0;
    double closestDvz_p16 = 9999;
    double closestDvz_p17 = 9999;
    for (int i = 0; i < p16id_event.vz.size(); ++i) {
      if (fabs(p16id_event.vz[i] - vz_p16) < closestDvz_p16 && i != id_p16) {
        closestDvz_p16 = fabs(p16id_event.vz[i] - vz_p16);
        vz_secondary_p16 = p16id_event.vz[i];
        nprim_secondary_p16 = p16id_event.nprim[i];
        dca_secondary_p16 = p16id_event.dca[i];
        pt_secondary_p16 = p16id_event.pt[i];
      }
    }
    for (int i = 0; i < p17id_event.vz.size(); ++i) {
      if (fabs(p17id_event.vz[i] - vz_p17) < closestDvz_p17 && i != id_p17) {
        closestDvz_p17 = fabs(p17id_event.vz[i] - vz_p17);
        vz_secondary_p17 = p17id_event.vz[i];
        nprim_secondary_p17 = p17id_event.nprim[i];
        dca_secondary_p17 = p17id_event.dca[i];
        pt_secondary_p17 = p17id_event.pt[i];
      }
    }
    
    dVzToNominal_p16->Fill(cent16, fabs(vz_p16 - vz_secondary_p16));
    closestDCA_p16->Fill(cent16, dca_secondary_p16);
    closestPt_p16->Fill(cent16, pt_secondary_p16);
    closestMult_p16->Fill(cent16, nprim_secondary_p16);
    
    dVzToNominal_p17->Fill(cent17, fabs(vz_p17 - vz_secondary_p17));
    closestDCA_p17->Fill(cent17, dca_secondary_p17);
    closestPt_p17->Fill(cent17, pt_secondary_p17);
    closestMult_p17->Fill(cent17, nprim_secondary_p17);
    
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
  
  rank_p16id->Scale(1.0 / rank_p16id->Integral());
  rank_p17id->Scale(1.0 / rank_p17id->Integral());
  Vz_p16id->Scale(1.0 / Vz_p16id->Integral());
  Vz_p17id->Scale(1.0 / Vz_p17id->Integral());
  dVz_p16id->Scale(1.0 / dVz_p16id->Integral());
  dVz_p17id->Scale(1.0 / dVz_p17id->Integral());
  dVz->Scale(1.0 / dVz->Integral());
  dnprim->Scale(1.0 / dnprim->Integral());
  dnglobal->Scale(1.0 / dnglobal->Integral());
  dnprimmatched->Scale(1.0 / dnprimmatched->Integral());
  dVzmatched->Scale(1.0 / dVzmatched->Integral());

  Overlay1D(rank_p16id, rank_p17id, "P16id", "P17id", hopts, copts, FLAGS_outdir, "compare_rank", "", "rank", "fraction");

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
 
  TProfile * nomP16Mult_avg = (TProfile*) nominalMult_p16->ProfileX();
  TProfile * secP16Mult_avg = (TProfile*) closestMult_p16->ProfileX();
  TProfile * nomP17Mult_avg = (TProfile*) nominalMult_p17->ProfileX();
  TProfile * secP17Mult_avg = (TProfile*) closestMult_p17->ProfileX();

  TProfile * nomP16DCA_avg = (TProfile*) nominalDCA_p16->ProfileX();
  TProfile * secP16DCA_avg = (TProfile*) closestDCA_p16->ProfileX();
  TProfile * nomP17DCA_avg = (TProfile*) nominalDCA_p17->ProfileX();
  TProfile * secP17DCA_avg = (TProfile*) closestDCA_p17->ProfileX();
  
  TProfile * nomP16Pt_avg = (TProfile*) nominalPt_p16->ProfileX();
  TProfile * secP16Pt_avg = (TProfile*) closestPt_p16->ProfileX();
  TProfile * nomP17Pt_avg = (TProfile*) nominalPt_p17->ProfileX();
  TProfile * secP17Pt_avg = (TProfile*) closestPt_p17->ProfileX();
  
  TProfile * p16dvz_avg = (TProfile*) dVzToNominal_p16->ProfileX();
  TProfile * p17dvz_avg = (TProfile*) dVzToNominal_p17->ProfileX();
  
  p16dvz_avg->GetYaxis()->SetRangeUser(6, 24);
  Overlay1D( p16dvz_avg, p17dvz_avg, "P16id", "P17id", hopts, copts, FLAGS_outdir, "dvz_avg", "",
            "centrality", "<dV_{z} [cm]>");
  Overlay1D( nomP16Mult_avg, nomP17Mult_avg, "P16id nominal", "P17id nominal", hopts, copts, FLAGS_outdir, "primary_mult_avg", "",
            "centrality", "<N_{primary}>");
  secP16Mult_avg->GetYaxis()->SetRangeUser(0, 100);
  Overlay1D( secP16Mult_avg, secP17Mult_avg, "P16id secondary", "P17id secondary", hopts, copts, FLAGS_outdir, "secondary_mult_avg", "",
            "centrality", "<N_{primary}>");
  nomP16DCA_avg->GetYaxis()->SetRangeUser(0, 5);
  Overlay1D( nomP16DCA_avg, nomP17DCA_avg, "P16id nominal", "P17id nominal", hopts, copts, FLAGS_outdir, "primary_dca_avg", "",
            "centrality", "<DCA [cm]>");
  secP16DCA_avg->GetYaxis()->SetRangeUser(0, 5);
  Overlay1D( secP16DCA_avg, secP17DCA_avg, "P16id secondary", "P17id secondary", hopts, copts, FLAGS_outdir, "secondary_dca_avg", "",
            "centrality", "<DCA [cm]>");
  nomP16Pt_avg->GetYaxis()->SetRangeUser(0, 2);
  Overlay1D( nomP16Pt_avg, nomP17Pt_avg, "P16id nominal", "P17id nominal", hopts, copts, FLAGS_outdir, "primary_pt_avg", "",
            "centrality", "<p_{T}>");
  secP16Pt_avg->GetYaxis()->SetRangeUser(0, 2);
  Overlay1D( secP16Pt_avg, secP17Pt_avg, "P16id secondary", "P17id secondary", hopts, copts, FLAGS_outdir, "secondary_pt_avg", "",
            "centrality", "<p_{T}>");
  
  TH1D* central_dvz_p16 = (TH1D*) dVzToNominal_p16->ProjectionY("tmp111", 1, 1);
  TH1D* central_dvz_p17 = (TH1D*) dVzToNominal_p17->ProjectionY("tmp222", 1, 1);
  
  central_dvz_p16->Scale(1.0 / central_dvz_p16->Integral());
  central_dvz_p17->Scale(1.0 / central_dvz_p17->Integral());
  central_dvz_p16->GetYaxis()->SetRangeUser(0.0, 0.2);
  Overlay1D(central_dvz_p16, central_dvz_p17, "P16id", "P17id", hopts, copts, FLAGS_outdir, "dvz_central", "", "dV_{z}", "fraction");
  
  gflags::ShutDownCommandLineFlags();
  return 0;
}

//template<typename H>
// void Print2DSimple(H* h,
//                    histogramOpts hopts,
//                    canvasOpts copts,
//                    std::string output_loc,
//                    std::string output_name,
//                    std::string canvas_title,
//                    std::string x_axis_label,
//                    std::string y_axis_label,
//                    std::string opt = "COLZ")

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
