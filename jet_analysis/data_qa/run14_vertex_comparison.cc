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
#include "TProfile2D.h"
#include "TF1.h"
#include "TStyle.h"

#include "TCanvas.h"

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

// compare P16id, P17id, P18if on an event-by-event level
// from trees produced by the muDsts

DEFINE_string(name, "compare_y14", "output name");
DEFINE_string(d1, "p16id.root", "input root file for P16id");
DEFINE_string(d1Name, "p16id", "data source 1 string identifier");
DEFINE_string(d2, "p17id.root", "input root file for P17id");
DEFINE_string(d2Name, "p17id", "data source 2 string identifier");
DEFINE_string(d3, "p18if.root", "input root file for P18if");
DEFINE_string(d3Name, "p18if", "data source 3 string identifier");
DEFINE_string(tree, "tree", "input tree name");
DEFINE_string(outdir, "compare_y14", "output directory");
DEFINE_double(dVzMax, 3.0, "maximum vz - vpdVz");

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
  std::vector<std::vector<double>> pt;
  std::vector<std::vector<double>> dca;
  std::vector<std::vector<double>> eta;
  std::vector<std::vector<double>> nhit;
  std::vector<std::vector<double>> nhitposs;
  std::vector<std::vector<double>> hft;
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
    TTreeReaderValue<std::vector<std::vector<double>>> pt(reader, "pt");
    TTreeReaderValue<std::vector<std::vector<double>>> dca(reader, "dca");
    TTreeReaderValue<std::vector<std::vector<double>>> eta(reader, "eta");
    TTreeReaderValue<std::vector<std::vector<double>>> nhit(reader, "nhit");
    TTreeReaderValue<std::vector<std::vector<double>>> nhitposs(reader, "nhitposs");
    TTreeReaderValue<std::vector<std::vector<double>>> hft(reader, "hft");
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
      evt.eta = *eta;
      evt.nhit = *nhit;
      evt.nhitposs = *nhitposs;
      evt.hft = *hft;
      map[{*runid, *eventid}] = evt;
    }
    
    return true;
  }
  else
    return false;
}

int SelectVertex(const Event& event, double dVzRange) {
  if (event.vz.size() == 0)
    return -1;
  double vpdvz = event.vpdvz;
  for (int i = 0; i < event.vz.size(); ++i)
    if (fabs(vpdvz - event.vz[i]) < dVzRange)
      return i;
  return -1;
}

int SelectSecondaryVertex(const Event& event, int primary_idx) {
  if (event.vz.size() == 0)
    return -1;
  int maxNPrim = -1;
  int secondaryIdx = -1;
  for (int i = 0; i < event.vz.size(); ++i)
    if (i != primary_idx && event.nprim[i] > maxNPrim) {
      secondaryIdx = i;
      maxNPrim = event.nprim[i];
    }
  
  return secondaryIdx;
}

int EventLoop(const Event& map, int vtx/*, TH1D* h_nprim*/, TH1D* h_pt) {

  const vector<double>& pt = map.pt.at(vtx);
  const vector<double>& eta = map.eta.at(vtx);
  const vector<double>& dca = map.dca.at(vtx);
  const vector<double>& nhit = map.nhit.at(vtx);
  const vector<double>& nhitposs = map.nhitposs.at(vtx);
  const vector<double>& hft = map.hft.at(vtx);

  int nprim = 0;
  for (int i = 0; i < pt.size(); ++i) {

    if (hft.at(i) > 0)
      LOG(INFO) << "WTF";

    if (fabs(eta.at(i)) > 1.0)
      continue;
    if (nhit.at(i) < 20)
      continue;
    if (dca.at(i) > 3.0)
      continue;

    nprim++;
    h_pt->Fill(pt.at(i));
  }
  
  return nprim;
}

void FillInclusive(Event& event, int vtx_id, TH1D* vz, TH1D* dvz, TH1D* nprim) {
  for(int i = 0; i < event.vz.size(); ++i) {
    if (i == vtx_id)
      continue;
    vz->Fill(event.vz.at(i));
    dvz->Fill(event.vz.at(vtx_id) - event.vz.at(i));
    nprim->Fill(event.nprim.at(i));
  }
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
  
  // define an enumeration for our sources, and hard-code
  // the number of
  enum SOURCE {
    d1 = 0,
    d2 = 1,
    d3 = 2
  };
  const vector<string> SOURCENAME {FLAGS_d1Name, FLAGS_d2Name, FLAGS_d3Name};
  
  // check to make sure the input file exists
  if (!boost::filesystem::exists(FLAGS_d1)) {
    std::cerr << "input file does not exist: " << FLAGS_d2 << std::endl;;
    return 1;
  }
  if (!boost::filesystem::exists(FLAGS_d2)) {
    std::cerr << "input file does not exist: " << FLAGS_d2 << std::endl;;
    return 1;
  }
  
  // build output directory if it doesn't exist, using boost::filesystem
  if (FLAGS_outdir.empty())
    FLAGS_outdir = "tmp";
  boost::filesystem::path dir(FLAGS_outdir.c_str());
  boost::filesystem::create_directories(dir);
  
  
  event_map d1_events;
  event_map d2_events;
  event_map d3_events;
  LOG(INFO) << "loading " << FLAGS_d1Name << " tree";
  ReadInTree(FLAGS_d1, FLAGS_tree, d1_events);
  LOG(INFO) << "loading " << FLAGS_d2Name << " tree";
  ReadInTree(FLAGS_d2, FLAGS_tree, d2_events);
  LOG(INFO) << "loading " << FLAGS_d3Name << " tree";
  ReadInTree(FLAGS_d3, FLAGS_tree, d3_events);
  LOG(INFO) << "successfully loaded in trees";
 
  // histograms
  // ----------
  
  // vertex rank
  vector<TH1D*> h_nvertices;
  vector<TH1D*> h_rank;
  vector<TH1D*> h_vz;
  vector<TH1D*> h_vpdvz;
  vector<TH1D*> h_dvz;
  vector<TH1D*> h_nprim;
  vector<TH1D*> h_nprim_ana;
  vector<TH1D*> h_nglob;
  vector<TH1D*> h_pt;
  vector<TH1D*> h_secondary_vz;
  vector<TH1D*> h_secondary_dvz;
  vector<TH1D*> h_secondary_nprim;
  vector<TH1D*> h_secondary_dn;
  
  vector<TH1D*> h_inclusive_vz;
  vector<TH1D*> h_inclusive_dvz;
  vector<TH1D*> h_inclusive_nprim;
  
  for (auto& prefix : SOURCENAME) {
    h_nvertices.push_back(new TH1D(MakeString(prefix, "nvert").c_str(), "", 50, 0, 50));
    h_rank.push_back(new TH1D(MakeString(prefix, "rank").c_str(), "", 50, -10, 3));
    h_vz.push_back(new TH1D(MakeString(prefix, "vz").c_str(), "", 70, -35, 35));
    h_vpdvz.push_back(new TH1D(MakeString(prefix, "vpdvz").c_str(), "", 70, -35, 35));
    h_dvz.push_back(new TH1D(MakeString(prefix, "dvz").c_str(), "", 50, -5, 5));
    h_nprim.push_back(new TH1D(MakeString(prefix, "nprim").c_str(), "", 100, 0, 1500));
    h_nprim_ana.push_back(new TH1D(MakeString(prefix, "nprimana").c_str(), "", 100, 0, 1000));
    h_nglob.push_back(new TH1D(MakeString(prefix, "nglob").c_str(), "", 100, 0, 5000));
    h_pt.push_back(new TH1D(MakeString(prefix, "pt").c_str(), "", 100, 0, 5));
    h_secondary_vz.push_back(new TH1D(MakeString(prefix, "secondaryvz").c_str(), "", 100, -100, 100));
    h_secondary_dvz.push_back(new TH1D(MakeString(prefix, "secondarydvz").c_str(), "", 100, -50, 50));
    h_secondary_nprim.push_back(new TH1D(MakeString(prefix, "secondarynprim").c_str(), "", 100, 0, 1500));
    h_secondary_dn.push_back(new TH1D(MakeString(prefix, "secondarydn").c_str(), "", 200, -1500, 1500));
    h_inclusive_vz.push_back(new TH1D(MakeString(prefix, "inclusivevz").c_str(), "", 100, -100, 100));
    h_inclusive_dvz.push_back(new TH1D(MakeString(prefix, "inclusivedvz").c_str(), "", 100, -50, 50));
    h_inclusive_nprim.push_back(new TH1D(MakeString(prefix, "inclusivenprim").c_str(), "", 100, 0, 1500));
  }
  
  TH1D* dvz_d1_d2 = new TH1D("d1d2dvz", "", 100, -3, 3);
  TH1D* dvz_d1_d3 = new TH1D("d1d3dvz", "", 100, -3, 3);
  TH1D* dvz_d2_d3 = new TH1D("d2d3dvz", "", 100, -3, 3);
  
  TH1D* dn_d1_d2_same = new TH1D("d1d2dnsame", "", 100, 0, 1000);
  TH1D* dn_d1_d2_diff = new TH1D("d1d2dndiff", "", 100, 0, 1000);
  TH1D* dglob_d1_d2_same = new TH1D("d1d2dglobsame", "", 100, 0, 1000);
  TH1D* dglob_d1_d2_diff = new TH1D("d1d2dglobdiff", "", 100, 0, 1000);
  
  TH1D* dn_d1_d3_same = new TH1D("d1d3dnsame", "", 100, 0, 1000);
  TH1D* dn_d1_d3_diff = new TH1D("d1d3dndiff", "", 100, 0, 1000);
  TH1D* dglob_d1_d3_same = new TH1D("d1d3dglobsame", "", 100, 0, 1000);
  TH1D* dglob_d1_d3_diff = new TH1D("d1d3dglobdiff", "", 100, 0, 1000);
  
  TH1D* dn_d2_d3_same = new TH1D("d2d3dnsame", "", 100, 0, 1000);
  TH1D* dn_d2_d3_diff = new TH1D("d2d3dndiff", "", 100, 0, 1000);
  TH1D* dglob_d2_d3_same = new TH1D("d2d3dglobsame", "", 100, 0, 1000);
  TH1D* dglob_d2_d3_diff = new TH1D("d2d3dglobdiff", "", 100, 0, 1000);
  
  TProfile2D* dvz_d1_d2_mult = new TProfile2D("d1d2dvzmult", ";dV_{z};n_{primary}", 100, -6, 6, 100, 0, 1000);
  TProfile2D* dvz_d1_d3_mult = new TProfile2D("d1d3dvzmult", ";dV_{z};n_{primary}", 100, -6, 6, 100, 0, 1000);
  TProfile2D* dvz_d2_d3_mult = new TProfile2D("d2d3dvzmult", ";dV_{z};n_{primary}", 100, -6, 6, 100, 0, 1000);
  
  TProfile* mult_dvz_d1_d2 = new TProfile("multdvzd1d2", ";N_{primary};dV_{z}", 50, 0, 100);
  TProfile* mult_dvz_d1_d3 = new TProfile("multdvzd1d3", ";N_{primary};dV_{z}", 50, 0, 100);
  TProfile* mult_dvz_d2_d3 = new TProfile("multdvzd2d3", ";N_{primary};dV_{z}", 50, 0, 100);
  
  TProfile* mult_dmult_d1_d2 = new TProfile("multdmultd1d2", ";N_{primary};dN_{primary}", 50, 0, 100);
  TProfile* mult_dmult_d1_d3 = new TProfile("multdmultd1d3", ";N_{primary};dN_{primary}", 50, 0, 100);
  TProfile* mult_dmult_d2_d3 = new TProfile("multdmultd2d3", ";N_{primary};dN_{primary}", 50, 0, 100);
  
  for (auto& entry : d1_events) {
  
    auto key = entry.first;
    auto d1_event = entry.second;
    
    // make sure there is no hft
    if (d1_event.pxl != 0 || d1_event.ist != 0)
      continue;
    
    // check if the event is present in the other datasets
    if (d2_events.find(key) == d2_events.end() ||
        d3_events.find(key) == d3_events.end())
      continue;
    auto d2_event = d2_events[key];
    auto d3_event = d3_events[key];
    
    // check there is no hft in the other datasets
    if (d2_event.pxl != 0 || d2_event.ist != 0 ||
        d3_event.pxl != 0 || d3_event.ist != 0)
      continue;
    
    // select vertex for all three
    int d1_vtx_id = SelectVertex(d1_event, FLAGS_dVzMax);
    int d2_vtx_id = SelectVertex(d2_event, FLAGS_dVzMax);
    int d3_vtx_id = SelectVertex(d3_event, FLAGS_dVzMax);
    
    if (d1_vtx_id < 0 || d2_vtx_id < 0 || d3_vtx_id < 0)
      continue;
    
    // plot inclusive (minus primary vertex) quantities
    FillInclusive(d1_event, d1_vtx_id, h_inclusive_vz[0], h_inclusive_dvz[0], h_inclusive_nprim[0]);
    FillInclusive(d2_event, d2_vtx_id, h_inclusive_vz[1], h_inclusive_dvz[1], h_inclusive_nprim[1]);
    FillInclusive(d3_event, d3_vtx_id, h_inclusive_vz[2], h_inclusive_dvz[2], h_inclusive_nprim[2]);
    
    // select secondary vertex for all three
    int d1_vtx_secondary = SelectSecondaryVertex(d1_event, d1_vtx_id);
    int d2_vtx_secondary = SelectSecondaryVertex(d2_event, d2_vtx_id);
    int d3_vtx_secondary = SelectSecondaryVertex(d3_event, d3_vtx_id);
    
    // plot secondary quanities
    if (d1_vtx_secondary >= 0) {
      h_secondary_vz[0]->Fill(d1_event.vz[d1_vtx_secondary]);
      h_secondary_dvz[0]->Fill(d1_event.vz[d1_vtx_id] - d1_event.vz[d1_vtx_secondary]);
      h_secondary_nprim[0]->Fill(d1_event.nprim[d1_vtx_secondary]);
      h_secondary_dn[0]->Fill(d1_event.nprim[d1_vtx_id] - d1_event.nprim[d1_vtx_secondary]);
    }
    
    if (d2_vtx_secondary >= 0) {
      h_secondary_vz[1]->Fill(d2_event.vz[d2_vtx_secondary]);
      h_secondary_dvz[1]->Fill(d2_event.vz[d2_vtx_id] - d2_event.vz[d2_vtx_secondary]);
      h_secondary_nprim[1]->Fill(d2_event.nprim[d2_vtx_secondary]);
      h_secondary_dn[1]->Fill(d2_event.nprim[d2_vtx_id] - d2_event.nprim[d2_vtx_secondary]);
    }
    
    if (d3_vtx_secondary >= 0) {
      h_secondary_vz[2]->Fill(d3_event.vz[d3_vtx_secondary]);
      h_secondary_dvz[2]->Fill(d3_event.vz[d3_vtx_id] - d3_event.vz[d3_vtx_secondary]);
      h_secondary_nprim[2]->Fill(d3_event.nprim[d3_vtx_secondary]);
      h_secondary_dn[2]->Fill(d3_event.nprim[d3_vtx_id] - d3_event.nprim[d3_vtx_secondary]);
    }

    // calculate nprim_ana and pT spectra
    int nprim_d1 = EventLoop(d1_event, d1_vtx_id, h_pt[0]);
    int nprim_d2 = EventLoop(d2_event, d2_vtx_id, h_pt[1]);
    int nprim_d3 = EventLoop(d3_event, d3_vtx_id, h_pt[2]);
    h_nprim_ana[0]->Fill(nprim_d1);
    h_nprim_ana[1]->Fill(nprim_d2);
    h_nprim_ana[2]->Fill(nprim_d3);
    
    // fill event level
    
    h_nvertices[d1]->Fill(d1_event.rank.size());
    h_nvertices[d2]->Fill(d2_event.rank.size());
    h_nvertices[d3]->Fill(d3_event.rank.size());
    
    h_rank[d1]->Fill(d1_event.rank[d1_vtx_id]);
    h_rank[d2]->Fill(d2_event.rank[d2_vtx_id]);
    h_rank[d3]->Fill(d3_event.rank[d3_vtx_id]);
    
    h_vz[d1]->Fill(d1_event.vz[d1_vtx_id]);
    h_vz[d2]->Fill(d2_event.vz[d2_vtx_id]);
    h_vz[d3]->Fill(d3_event.vz[d3_vtx_id]);
    
    h_vpdvz[d1]->Fill(d1_event.vpdvz);
    h_vpdvz[d2]->Fill(d2_event.vpdvz);
    h_vpdvz[d3]->Fill(d3_event.vpdvz);
    
    h_dvz[d1]->Fill(d1_event.vpdvz - d1_event.vz[d1_vtx_id]);
    h_dvz[d2]->Fill(d2_event.vpdvz - d2_event.vz[d2_vtx_id]);
    h_dvz[d3]->Fill(d3_event.vpdvz - d3_event.vz[d3_vtx_id]);
    
    h_nprim[d1]->Fill(d1_event.nprim[d1_vtx_id]);
    h_nprim[d2]->Fill(d2_event.nprim[d2_vtx_id]);
    h_nprim[d3]->Fill(d3_event.nprim[d3_vtx_id]);
    
    h_nglob[d1]->Fill(d1_event.nglobal);
    h_nglob[d2]->Fill(d2_event.nglobal);
    h_nglob[d3]->Fill(d3_event.nglobal);
    
    dvz_d1_d2->Fill(d1_event.vz[d1_vtx_id] - d2_event.vz[d2_vtx_id]);
    dvz_d1_d3->Fill(d1_event.vz[d1_vtx_id] - d3_event.vz[d3_vtx_id]);
    dvz_d2_d3->Fill(d2_event.vz[d2_vtx_id] - d3_event.vz[d3_vtx_id]);
    
    dvz_d1_d2_mult->Fill(d1_event.vz[d1_vtx_id] - d2_event.vz[d2_vtx_id], d1_event.nprim[d1_vtx_id], d1_event.nprim[d1_vtx_id] - d2_event.nprim[d2_vtx_id]);
    dvz_d1_d3_mult->Fill(d1_event.vz[d1_vtx_id] - d3_event.vz[d3_vtx_id], d1_event.nprim[d1_vtx_id], d1_event.nprim[d1_vtx_id] - d3_event.nprim[d3_vtx_id]);
    dvz_d2_d3_mult->Fill(d2_event.vz[d2_vtx_id] - d3_event.vz[d3_vtx_id], d1_event.nprim[d1_vtx_id], d2_event.nprim[d2_vtx_id] - d3_event.nprim[d3_vtx_id]);
    
    double d1_d2_dvz = fabs(d1_event.vz[d1_vtx_id] - d2_event.vz[d2_vtx_id]);
    int d1_d2_dnprim = fabs(d1_event.nprim[d1_vtx_id] - d2_event.nprim[d2_vtx_id]);
    int d1_d2_dnglobal = fabs(d1_event.nglobal - d2_event.nglobal);
    double d1_d3_dvz = fabs(d1_event.vz[d1_vtx_id] - d3_event.vz[d3_vtx_id]);
    int d1_d3_dnprim = fabs(d1_event.nprim[d1_vtx_id] - d3_event.nprim[d3_vtx_id]);
    int d1_d3_dnglobal = fabs(d1_event.nglobal - d3_event.nglobal);
    double d2_d3_dvz = fabs(d2_event.vz[d2_vtx_id] - d3_event.vz[d3_vtx_id]);
    int d2_d3_dnprim = fabs(d2_event.nprim[d2_vtx_id] - d3_event.nprim[d3_vtx_id]);
    int d2_d3_dnglobal = fabs(d2_event.nglobal - d3_event.nglobal);
    
    // fill the differences
    mult_dvz_d1_d2->Fill(d1_event.nprim[d1_vtx_id], fabs(d1_d2_dvz));
    mult_dvz_d1_d3->Fill(d1_event.nprim[d1_vtx_id], fabs(d1_d3_dvz));
    mult_dvz_d2_d3->Fill(d2_event.nprim[d2_vtx_id], fabs(d2_d3_dvz));
    
    if (d1_event.nprim[d1_vtx_id] > 0)
      mult_dmult_d1_d2->Fill(d1_event.nprim[d1_vtx_id], (double)d1_d2_dnprim/d1_event.nprim[d1_vtx_id]);
    else
      LOG(INFO) << "d1 zero primaries";
    if (d1_event.nprim[d1_vtx_id] > 0)
      mult_dmult_d1_d3->Fill(d1_event.nprim[d1_vtx_id], (double)d1_d3_dnprim/d1_event.nprim[d1_vtx_id]);
    else
      LOG(INFO) << "d1 zero primaries";
    if (d2_event.nprim[d2_vtx_id] > 0)
      mult_dmult_d2_d3->Fill(d3_event.nprim[d3_vtx_id], (double)d2_d3_dnprim/d3_event.nprim[d3_vtx_id]);
    else
      LOG(INFO) << "d2 zero primaries";
    if (d1_d2_dvz > 1.0) {
      dn_d1_d2_diff->Fill(d1_d2_dnprim);
      dglob_d1_d2_diff->Fill(d1_d2_dnglobal);
    }
    else {
      dn_d1_d2_same->Fill(d1_d2_dnprim);
      dglob_d1_d2_same->Fill(d1_d2_dnglobal);
    }
    if (d1_d3_dvz > 1.0) {
      dn_d1_d3_diff->Fill(d1_d3_dnprim);
      dglob_d1_d3_diff->Fill(d1_d3_dnglobal);
    }
    else {
      dn_d1_d3_same->Fill(d1_d3_dnprim);
      dglob_d1_d3_same->Fill(d1_d3_dnglobal);
    }
    if (d2_d3_dvz > 1.0) {
      dn_d2_d3_diff->Fill(d2_d3_dnprim);
      dglob_d2_d3_diff->Fill(d2_d3_dnglobal);
    }
    else {
      dn_d2_d3_same->Fill(d2_d3_dnprim);
      dglob_d2_d3_same->Fill(d2_d3_dnglobal);
    }
  }
  
  // scale
  for (int i = 0; i < SOURCENAME.size(); ++i) {
    h_nvertices[i]->Scale(1.0 / h_nvertices[i]->Integral());
    h_rank[i]->Scale(1.0 / h_rank[i]->Integral());
    h_vz[i]->Scale(1.0 / h_vz[i]->Integral());
    h_vpdvz[i]->Scale(1.0 / h_vpdvz[i]->Integral());
    h_dvz[i]->Scale(1.0 / h_dvz[i]->Integral());
    h_nprim[i]->Scale(1.0 / h_nprim[i]->Integral());
    h_nprim_ana[i]->Scale(1.0 / h_nprim_ana[i]->Integral());
    h_nglob[i]->Scale(1.0 / h_nglob[i]->Integral());
    h_pt[i]->Scale(1.0 / h_vz[i]->GetEntries());
    h_secondary_vz[i]->Scale(1.0 / h_secondary_vz[i]->Integral());
    h_secondary_dvz[i]->Scale(1.0 / h_secondary_dvz[i]->Integral());
    h_secondary_nprim[i]->Scale(1.0 / h_secondary_nprim[i]->Integral());
    h_secondary_dn[i]->Scale(1.0 / h_secondary_dn[i]->Integral());
    h_inclusive_vz[i]->Scale(1.0 / h_inclusive_vz[i]->Integral());
    h_inclusive_dvz[i]->Scale(1.0 / h_inclusive_dvz[i]->Integral());
    h_inclusive_nprim[i]->Scale(1.0 / h_inclusive_nprim[i]->Integral());
    
    dn_d1_d2_same->Scale(1.0 / dn_d1_d2_same->Integral());
    dn_d1_d3_same->Scale(1.0 / dn_d1_d3_same->Integral());
    dn_d2_d3_same->Scale(1.0 / dn_d2_d3_same->Integral());
    dn_d1_d2_diff->Scale(1.0 / dn_d1_d2_diff->Integral());
    dn_d1_d3_diff->Scale(1.0 / dn_d1_d3_diff->Integral());
    dn_d2_d3_diff->Scale(1.0 / dn_d2_d3_diff->Integral());
    
    dglob_d1_d2_same->Scale(1.0 / dglob_d1_d2_same->Integral());
    dglob_d1_d3_same->Scale(1.0 / dglob_d1_d3_same->Integral());
    dglob_d2_d3_same->Scale(1.0 / dglob_d2_d3_same->Integral());
    dglob_d1_d2_diff->Scale(1.0 / dglob_d1_d2_diff->Integral());
    dglob_d1_d3_diff->Scale(1.0 / dglob_d1_d3_diff->Integral());
    dglob_d2_d3_diff->Scale(1.0 / dglob_d2_d3_diff->Integral());
  }
  
  // print results
  // create our histogram and canvas options
  histogramOpts hopts;
  canvasOpts copts;
  canvasOpts coptslogz;
  coptslogz.log_z = true;
  canvasOpts coptslogy;
  coptslogy.log_y = true;
  canvasOpts coptslogynoleg;
  coptslogynoleg.log_y = true;
  coptslogynoleg.do_legend = false;
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
  
  Overlay1D(h_nvertices, SOURCENAME, hopts, copts, FLAGS_outdir, "nvertices",
            "", "N_{vertices}", "fraction", "");
  Overlay1D(h_rank, SOURCENAME, hopts, cOptsTopLeftLeg, FLAGS_outdir, "rank",
            "", "rank", "fraction", "");
  Overlay1D(h_vz, SOURCENAME, hopts, copts, FLAGS_outdir, "vz",
            "", "V_{z}", "fraction", "");
  Overlay1D(h_vpdvz, SOURCENAME, hopts, copts, FLAGS_outdir, "vpdvz",
            "", "VPD V_{z}", "fraction", "");
  Overlay1D(h_dvz, SOURCENAME, hopts, coptslogy, FLAGS_outdir, "dvz",
            "", "dV_{z}", "fraction", "");
  Overlay1D(h_nprim, SOURCENAME, hopts, coptslogy, FLAGS_outdir, "nprim",
            "", "N_{primary}", "fraction", "");
  Overlay1D(h_nprim_ana, SOURCENAME, hopts, coptslogy, FLAGS_outdir, "nprim",
            "", "N_{primary}", "fraction", "");
  Overlay1D(h_nglob, SOURCENAME, hopts, cOptsBottomLeftLegLogy, FLAGS_outdir, "nglob",
            "", "N_{global}", "fraction", "");
  Overlay1D(h_pt, SOURCENAME, hopts, coptslogy, FLAGS_outdir, "pt",
            "", "p_{T}", "fraction", "");
  Overlay1D(h_secondary_vz, SOURCENAME, hopts, copts, FLAGS_outdir, "secondaryvz",
            "", "V_{z}", "fraction", "");
  Overlay1D(h_secondary_dvz, SOURCENAME, hopts, copts, FLAGS_outdir, "secondarydvz",
            "", "dV_{z}", "fraction", "");
  Overlay1D(h_secondary_nprim, SOURCENAME, hopts, coptslogy, FLAGS_outdir, "secondarynprim",
            "", "N_{primary}", "fraction", "");
  Overlay1D(h_secondary_dn, SOURCENAME, hopts, copts, FLAGS_outdir, "secondarydn",
            "", "dN_{primary}", "fraction", "");
  Overlay1D(h_inclusive_vz, SOURCENAME, hopts, copts, FLAGS_outdir, "inclusivevz",
            "", "V_{z}", "fraction", "");
  Overlay1D(h_inclusive_dvz, SOURCENAME, hopts, copts, FLAGS_outdir, "inclusivedvz",
            "", "dV_{z}", "fraction", "");
  Overlay1D(h_inclusive_nprim, SOURCENAME, hopts, coptslogy, FLAGS_outdir, "inclusivenprim",
            "", "N_{primary}", "fraction", "");
  
  
  PrettyPrint1D(dvz_d1_d2, hopts, coptslogynoleg, "", FLAGS_outdir, "dvzd1d2", "", "dV_{z}", "fraction");
  PrettyPrint1D(dvz_d1_d3, hopts, coptslogynoleg, "", FLAGS_outdir, "dvzd1d3", "", "dV_{z}", "fraction");
  PrettyPrint1D(dvz_d2_d3, hopts, coptslogynoleg, "", FLAGS_outdir, "dvzd2d3", "", "dV_{z}", "fraction");
  
  Print2DSimple(dvz_d1_d2_mult, hopts, copts, FLAGS_outdir, "dvzd1d2mult", "", "dV_{z}", "N_{primary}");
  Print2DSimple(dvz_d1_d3_mult, hopts, copts, FLAGS_outdir, "dvzd1d3mult", "", "dV_{z}", "N_{primary}");
  Print2DSimple(dvz_d2_d3_mult, hopts, copts, FLAGS_outdir, "dvzd2d3mult", "", "dV_{z}", "N_{primary}");
  
  TProfile* test1 = dvz_d1_d2_mult->ProfileX();
  TProfile* test2 = dvz_d1_d2_mult->ProfileY();
  
  TProfile* test3 = dvz_d1_d3_mult->ProfileX();
  TProfile* test4 = dvz_d1_d3_mult->ProfileY();
  
  Overlay1D(test1, test3, MakeString(FLAGS_d1Name, "-", FLAGS_d2Name), MakeString(FLAGS_d1Name, "-", FLAGS_d3Name), hopts, copts, FLAGS_outdir, "dvzdN", "", "dV_{z}", "<dN_{primary}>");
  Overlay1D(test2, test4, MakeString(FLAGS_d1Name, "-", FLAGS_d2Name), MakeString(FLAGS_d1Name, "-", FLAGS_d3Name), hopts, cOptsBottomLeg, FLAGS_outdir, "nprimdN", "", "N_{primary}", "<dN_{primary}>");
  
  
  LOG(INFO) << "probability to find a vertex w/ higher multiplicity than VPD vertex P16id: " << h_secondary_dn[0]->Integral(1, h_secondary_dn[0]->GetXaxis()->FindBin(0.0));
  LOG(INFO) << "probability to find a vertex w/ higher multiplicity than VPD vertex P16id: " << h_secondary_dn[1]->Integral(1, h_secondary_dn[1]->GetXaxis()->FindBin(0.0));
  LOG(INFO) << "probability to find a vertex w/ higher multiplicity than VPD vertex P16id: " << h_secondary_dn[2]->Integral(1, h_secondary_dn[2]->GetXaxis()->FindBin(0.0));
  
  std::vector<string> diff_names{MakeString(SOURCENAME[0], "-", SOURCENAME[1]),
                                 MakeString(SOURCENAME[0], "-", SOURCENAME[2]),
                                 MakeString(SOURCENAME[1], "-", SOURCENAME[2])};
  std::vector<TH1D*> diff_nprim{dn_d1_d2_diff, dn_d1_d3_diff, dn_d2_d3_diff};
  std::vector<TH1D*> same_nprim{dn_d1_d2_same, dn_d1_d3_same, dn_d2_d3_same};
  std::vector<TH1D*> diff_nglob{dglob_d1_d2_diff, dglob_d1_d3_diff, dglob_d2_d3_diff};
  std::vector<TH1D*> same_nglob{dglob_d1_d2_same, dglob_d1_d3_same, dglob_d2_d3_same};
  
  Overlay1D(diff_nprim, diff_names, hopts, coptslogy, FLAGS_outdir, "diffdnprim", "", "dN_{primary}", "fraction");
  Overlay1D(same_nprim, diff_names, hopts, coptslogy, FLAGS_outdir, "samednprim", "", "dN_{primary}", "fraction");
  Overlay1D(diff_nglob, diff_names, hopts, coptslogy, FLAGS_outdir, "diffdnglob", "", "dN_{global}", "fraction");
  Overlay1D(same_nglob, diff_names, hopts, coptslogy, FLAGS_outdir, "samednglob", "", "dN_{global}", "fraction");
  
  std::vector<TProfile*> mult_dvz;
  mult_dvz.push_back(mult_dvz_d1_d2);
  mult_dvz.push_back(mult_dvz_d1_d3);
  mult_dvz.push_back(mult_dvz_d2_d3);
  
  std::vector<TProfile*> mult_dmult;
  mult_dmult.push_back(mult_dmult_d1_d2);
  mult_dmult.push_back(mult_dmult_d1_d3);
  mult_dmult.push_back(mult_dmult_d2_d3);
  
  Overlay1D(mult_dvz, diff_names, hopts, copts, FLAGS_outdir, "multdvz", "", "N_{primary}", "|dV_{z}|", "", true);
  Overlay1D(mult_dmult, diff_names, hopts, copts, FLAGS_outdir, "multdmult", "", "N_{primary}", "|dN_{primary}|/N_{primary}", "", true);
  
  
//
//  void Print2DSimple(H* h,
//                     histogramOpts hopts,
//                     canvasOpts copts,
//                     std::string output_loc,
//                     std::string output_name,
//                     std::string canvas_title,
//                     std::string x_axis_label,
//                     std::string y_axis_label,
//                     std::string opt = "COLZ")
  
//  void PrettyPrint1D(H* h,
//                     histogramOpts hopts,
//                     canvasOpts copts,
//                     std::string hist_title,
//                     std::string output_loc,
//                     std::string output_name,
//                     std::string canvas_title,
//                     std::string x_axis_label,
//                     std::string y_axis_label,
//                     std::string legend_title = "") {
  
//  TH1D* rank_p16id = new TH1D("p16id_rank", "", 150, -10, 5);
//  TH1D* Vz_p16id = new TH1D("p16id_vz", "", 100, -50, 50);
//  TH1D* dVz_p16id = new TH1D("p16id_dvz", "", 100, 0, 50);
//  TH1D* Vz_p17id = new TH1D("p17id_vz", "", 100, -50, 50);
//  TH1D* rank_p17id = new TH1D("p17id_rank", "", 150, -10, 5);
//  TH1D* dVz_p17id = new TH1D("p17id_dvz", "", 100, 0, 50);
//  TH1D* dVz = new TH1D("dvz", "", 100, 0, 20);
//  TH1D* dVzmatched = new TH1D("dvzmatched", "", 100, 0, 20);
//  TH1D* dnprim = new TH1D("dnprim", "", 200, -500, 500);
//  TH1D* dnprimmatched = new TH1D("dnprimmatched", "", 200, -500, 500);
//  TH1D* dnglobal = new TH1D("dnglobal", "", 200, -500, 500);
//
//  // some centrality analyses
//  TH2D* dVzToNominal_p16 = new TH2D("dvztonom16", "", 8, 0, 8, 100, 0, 50);
//  TH2D* nominalDCA_p16   = new TH2D("nominaldca16", "", 8, 0, 8, 100, 0, 5);
//  TH2D* nominalPt_p16    = new TH2D("nominalpt16", "", 8, 0, 8, 100, 0, 2);
//  TH2D* nominalMult_p16  = new TH2D("nominalmult16", "", 8, 0, 8, 100, 0, 1500);
//  TH2D* closestDCA_p16   = new TH2D("closestdca16", "", 8, 0, 8, 100, 0, 5);
//  TH2D* closestPt_p16    = new TH2D("closestpt16", "", 8, 0, 8, 100, 0, 2);
//  TH2D* closestMult_p16  = new TH2D("closestmult16", "", 8, 0, 8, 100, 0, 1500);
//
//  TH2D* dVzToNominal_p17 = new TH2D("dvztonom17", "", 8, 0, 8, 100, 0, 50);
//  TH2D* nominalDCA_p17   = new TH2D("nominaldca17", "", 8, 0, 8, 100, 0, 5);
//  TH2D* nominalPt_p17    = new TH2D("nominalpt17", "", 8, 0, 8, 100, 0, 2);
//  TH2D* nominalMult_p17  = new TH2D("nominalmult17", "", 8, 0, 8, 100, 0, 1500);
//  TH2D* closestDCA_p17   = new TH2D("closestdca17", "", 8, 0, 8, 100, 0, 5);
//  TH2D* closestPt_p17    = new TH2D("closestpt17", "", 8, 0, 8, 100, 0, 2);
//  TH2D* closestMult_p17  = new TH2D("closestmult17", "", 8, 0, 8, 100, 0, 1500);
//
//  // create centrality definition
//  CentralityRun14 centrality;
//
//  for (auto& entry : p16id) {
//    auto key = entry.first;
//    auto p16id_event = entry.second;
//    if (p16id_event.pxl != 0 || p16id_event.ist != 0)
//      continue;
//    if (p17id.find(key) == p17id.end())
//      continue;
//    auto p17id_event = p17id[key];
//
//    int id_p16 = -1;
//    int nprim_p16 = -1;
//    double vz_p16 = 0.0;
//    double dca_p16 = 0;
//    double pt_p16 = 0;
//    int id_p17 = -1;
//    int nprim_p17 = -1;
//    double vz_p17 = 0.0;
//    double dca_p17 = 0;
//    double pt_p17 = 0;
//
//    for (int i = 0; i < p16id_event.vz.size(); ++i) {
//      if (fabs(p16id_event.vz[i] - p16id_event.vpdvz) < 3.0) {
//        id_p16 = i;
//        vz_p16 = p16id_event.vz[i];
//        nprim_p16 = p16id_event.nprim[i];
//        dca_p16 = p16id_event.dca[i];
//        pt_p16 = p16id_event.pt[i];
//        break;
//      }
//    }
//    for (int i = 0; i < p17id_event.vz.size(); ++i) {
//      if (fabs(p17id_event.vz[i] - p17id_event.vpdvz) < 3.0) {
//        id_p17 = i;
//        vz_p17 = p17id_event.vz[i];
//        nprim_p17 = p17id_event.nprim[i];
//        dca_p17 = p17id_event.dca[i];
//        pt_p17 = p17id_event.pt[i];
//        break;
//      }
//    }
//
//    if (nprim_p16 < 0 || nprim_p17 < 0)
//      continue;
//
//    rank_p16id->Fill(p16id_event.rank[id_p16]);
//    rank_p17id->Fill(p17id_event.rank[id_p17]);
//
//    Vz_p16id->Fill(p16id_event.vz[id_p16]);
//    Vz_p17id->Fill(p17id_event.vz[id_p17]);
//
//    dVz_p16id->Fill(fabs(p16id_event.vz[id_p16] - p16id_event.vpdvz));
//    dVz_p17id->Fill(fabs(p17id_event.vz[id_p17] - p17id_event.vpdvz));
//
//    dVz->Fill(fabs(p16id_event.vz[id_p16] - p17id_event.vz[id_p17]));
//    dnprim->Fill(p16id_event.nprim[id_p16] - p17id_event.vz[id_p17]);
//    dnglobal->Fill(p16id_event.nglobal - p17id_event.nglobal);
//
//
//
//    centrality.setEvent(p16id_event.runid, p16id_event.refmult[id_p16], p16id_event.zdcrate, vz_p16);
//    int cent16 = centrality.centrality9();
//    centrality.setEvent(p17id_event.runid, p17id_event.refmult[id_p17], p17id_event.zdcrate, vz_p17);
//    int cent17 = centrality.centrality9();
//
//    nominalDCA_p16->Fill(cent16, dca_p16);
//    nominalPt_p16->Fill(cent16, pt_p16);
//    nominalMult_p16->Fill(cent16, nprim_p16);
//
//    nominalDCA_p17->Fill(cent17, dca_p17);
//    nominalPt_p17->Fill(cent17, pt_p17);
//    nominalMult_p17->Fill(cent17, nprim_p17);
//
//    dnprimmatched->Fill(nprim_p16 - nprim_p17);
//    dVzmatched->Fill(vz_p16 - vz_p17);
//
//    // find the secondary vertices
//    int nprim_secondary_p16 = -1;
//    int nprim_secondary_p17 = -1;
//    double vz_secondary_p16 = 0;
//    double vz_secondary_p17 = 0;
//    double dca_secondary_p16 = 0;
//    double dca_secondary_p17 = 0;
//    double pt_secondary_p16 = 0;
//    double pt_secondary_p17 = 0;
//    double closestDvz_p16 = 9999;
//    double closestDvz_p17 = 9999;
//    for (int i = 0; i < p16id_event.vz.size(); ++i) {
//      if (fabs(p16id_event.vz[i] - vz_p16) < closestDvz_p16 && i != id_p16) {
//        closestDvz_p16 = fabs(p16id_event.vz[i] - vz_p16);
//        vz_secondary_p16 = p16id_event.vz[i];
//        nprim_secondary_p16 = p16id_event.nprim[i];
//        dca_secondary_p16 = p16id_event.dca[i];
//        pt_secondary_p16 = p16id_event.pt[i];
//      }
//    }
//    for (int i = 0; i < p17id_event.vz.size(); ++i) {
//      if (fabs(p17id_event.vz[i] - vz_p17) < closestDvz_p17 && i != id_p17) {
//        closestDvz_p17 = fabs(p17id_event.vz[i] - vz_p17);
//        vz_secondary_p17 = p17id_event.vz[i];
//        nprim_secondary_p17 = p17id_event.nprim[i];
//        dca_secondary_p17 = p17id_event.dca[i];
//        pt_secondary_p17 = p17id_event.pt[i];
//      }
//    }
//
//    dVzToNominal_p16->Fill(cent16, fabs(vz_p16 - vz_secondary_p16));
//    closestDCA_p16->Fill(cent16, dca_secondary_p16);
//    closestPt_p16->Fill(cent16, pt_secondary_p16);
//    closestMult_p16->Fill(cent16, nprim_secondary_p16);
//
//    dVzToNominal_p17->Fill(cent17, fabs(vz_p17 - vz_secondary_p17));
//    closestDCA_p17->Fill(cent17, dca_secondary_p17);
//    closestPt_p17->Fill(cent17, pt_secondary_p17);
//    closestMult_p17->Fill(cent17, nprim_secondary_p17);
//
//  }
//
//  // print results
//  // create our histogram and canvas options
//  histogramOpts hopts;
//  canvasOpts copts;
//  canvasOpts coptslogz;
//  coptslogz.log_z = true;
//  canvasOpts coptslogy;
//  coptslogy.log_y = true;
//  canvasOpts cOptsBottomLeg;
//  cOptsBottomLeg.leg_upper_bound = 0.4;
//  cOptsBottomLeg.leg_lower_bound = 0.18;
//  cOptsBottomLeg.leg_right_bound = 0.9;
//  cOptsBottomLeg.leg_left_bound = 0.7;
//  canvasOpts cOptsBottomLeftLeg;
//  cOptsBottomLeftLeg.leg_upper_bound = 0.4;
//  cOptsBottomLeftLeg.leg_lower_bound = 0.18;
//  cOptsBottomLeftLeg.leg_right_bound = 0.18;
//  cOptsBottomLeftLeg.leg_left_bound = 0.4;
//  canvasOpts cOptsBottomLeftLegLogy;
//  cOptsBottomLeftLegLogy.log_y = true;
//  cOptsBottomLeftLegLogy.leg_upper_bound = 0.4;
//  cOptsBottomLeftLegLogy.leg_lower_bound = 0.18;
//  cOptsBottomLeftLegLogy.leg_right_bound = 0.18;
//  cOptsBottomLeftLegLogy.leg_left_bound = 0.4;
//  canvasOpts cOptsTopLeftLeg;
//  cOptsTopLeftLeg.leg_right_bound = 0.18;
//  cOptsTopLeftLeg.leg_left_bound = 0.4;
//
//  rank_p16id->Scale(1.0 / rank_p16id->Integral());
//  rank_p17id->Scale(1.0 / rank_p17id->Integral());
//  Vz_p16id->Scale(1.0 / Vz_p16id->Integral());
//  Vz_p17id->Scale(1.0 / Vz_p17id->Integral());
//  dVz_p16id->Scale(1.0 / dVz_p16id->Integral());
//  dVz_p17id->Scale(1.0 / dVz_p17id->Integral());
//  dVz->Scale(1.0 / dVz->Integral());
//  dnprim->Scale(1.0 / dnprim->Integral());
//  dnglobal->Scale(1.0 / dnglobal->Integral());
//  dnprimmatched->Scale(1.0 / dnprimmatched->Integral());
//  dVzmatched->Scale(1.0 / dVzmatched->Integral());
//
//  Overlay1D(rank_p16id, rank_p17id, FLAGS_d1Name, FLAGS_d2Name, hopts, copts, FLAGS_outdir, "compare_rank", "", "rank", "fraction");
//
//  Overlay1D(Vz_p16id, Vz_p17id, FLAGS_d1Name, FLAGS_d2Name, hopts, copts, FLAGS_outdir, "compare_vz",
//            "", "V_{z}", "fraction");
//  Overlay1D(dVz_p16id, dVz_p17id, FLAGS_d1Name, FLAGS_d2Name, hopts, coptslogy, FLAGS_outdir, "compare_dvz",
//            "", "V_{z}", "fraction");
//  PrettyPrint1D(dVz, hopts, coptslogy, "dV_{z}", FLAGS_outdir, "dvz", "", "dV_{z}", "fraction");
//  PrettyPrint1D(dnprim, hopts, copts, "dN_{primary}", FLAGS_outdir, "dnprim", "", "dN_{primary}", "fraction");
//  PrettyPrint1D(dnprimmatched, hopts, coptslogy, "dN_{primary}", FLAGS_outdir, "dnprimmatched", "", "dN_{primary}", "fraction");
//  PrettyPrint1D(dnglobal, hopts, coptslogy, "dN_{global}", FLAGS_outdir, "dnglobal", "", "dN_{global}", "fraction");
//  PrettyPrint1D(dVzmatched, hopts, coptslogy, "dV_{z}", FLAGS_outdir, "dvzmatched", "", "dV_{z}", "fraction");
//
//  LOG(INFO) << "probability to have 10 or more: " << dnglobal->Integral(dnglobal->GetXaxis()->FindBin(10), 200);
//  LOG(INFO) << "probability to have 100 or more: " << dnglobal->Integral(dnglobal->GetXaxis()->FindBin(100), 200);
//
//  TProfile * nomP16Mult_avg = (TProfile*) nominalMult_p16->ProfileX();
//  TProfile * secP16Mult_avg = (TProfile*) closestMult_p16->ProfileX();
//  TProfile * nomP17Mult_avg = (TProfile*) nominalMult_p17->ProfileX();
//  TProfile * secP17Mult_avg = (TProfile*) closestMult_p17->ProfileX();
//
//  TProfile * nomP16DCA_avg = (TProfile*) nominalDCA_p16->ProfileX();
//  TProfile * secP16DCA_avg = (TProfile*) closestDCA_p16->ProfileX();
//  TProfile * nomP17DCA_avg = (TProfile*) nominalDCA_p17->ProfileX();
//  TProfile * secP17DCA_avg = (TProfile*) closestDCA_p17->ProfileX();
//
//  TProfile * nomP16Pt_avg = (TProfile*) nominalPt_p16->ProfileX();
//  TProfile * secP16Pt_avg = (TProfile*) closestPt_p16->ProfileX();
//  TProfile * nomP17Pt_avg = (TProfile*) nominalPt_p17->ProfileX();
//  TProfile * secP17Pt_avg = (TProfile*) closestPt_p17->ProfileX();
//
//  TProfile * p16dvz_avg = (TProfile*) dVzToNominal_p16->ProfileX();
//  TProfile * p17dvz_avg = (TProfile*) dVzToNominal_p17->ProfileX();
//
//  p16dvz_avg->GetYaxis()->SetRangeUser(6, 24);
//  Overlay1D( p16dvz_avg, p17dvz_avg, FLAGS_d1Name, FLAGS_d2Name, hopts, copts, FLAGS_outdir, "dvz_avg", "",
//            "centrality", "<dV_{z} [cm]>");
//  Overlay1D( nomP16Mult_avg, nomP17Mult_avg, FLAGS_d1Name + " nominal", FLAGS_d2Name + " nominal", hopts, copts, FLAGS_outdir, "primary_mult_avg", "",
//            "centrality", "<N_{primary}>");
//  secP16Mult_avg->GetYaxis()->SetRangeUser(0, 100);
//  Overlay1D( secP16Mult_avg, secP17Mult_avg, FLAGS_d1Name + " secondary", FLAGS_d2Name + " secondary", hopts, copts, FLAGS_outdir, "secondary_mult_avg", "",
//            "centrality", "<N_{primary}>");
//  nomP16DCA_avg->GetYaxis()->SetRangeUser(0, 5);
//  Overlay1D( nomP16DCA_avg, nomP17DCA_avg, FLAGS_d1Name + " nominal", FLAGS_d2Name + " nominal", hopts, copts, FLAGS_outdir, "primary_dca_avg", "",
//            "centrality", "<DCA [cm]>");
//  secP16DCA_avg->GetYaxis()->SetRangeUser(0, 5);
//  Overlay1D( secP16DCA_avg, secP17DCA_avg, FLAGS_d1Name + " secondary", FLAGS_d2Name + " secondary", hopts, copts, FLAGS_outdir, "secondary_dca_avg", "",
//            "centrality", "<DCA [cm]>");
//  nomP16Pt_avg->GetYaxis()->SetRangeUser(0, 2);
//  Overlay1D( nomP16Pt_avg, nomP17Pt_avg, FLAGS_d1Name + " nominal", FLAGS_d2Name + " nominal", hopts, copts, FLAGS_outdir, "primary_pt_avg", "",
//            "centrality", "<p_{T}>");
//  secP16Pt_avg->GetYaxis()->SetRangeUser(0, 2);
//  Overlay1D( secP16Pt_avg, secP17Pt_avg, FLAGS_d1Name + " secondary", FLAGS_d2Name + " secondary", hopts, copts, FLAGS_outdir, "secondary_pt_avg", "",
//            "centrality", "<p_{T}>");
//
//  TH1D* central_dvz_p16 = (TH1D*) dVzToNominal_p16->ProjectionY("tmp111", 1, 1);
//  TH1D* central_dvz_p17 = (TH1D*) dVzToNominal_p17->ProjectionY("tmp222", 1, 1);
//
//  central_dvz_p16->Scale(1.0 / central_dvz_p16->Integral());
//  central_dvz_p17->Scale(1.0 / central_dvz_p17->Integral());
//  central_dvz_p16->GetYaxis()->SetRangeUser(0.0, 0.2);
//  Overlay1D(central_dvz_p16, central_dvz_p17, FLAGS_d1Name, FLAGS_d2Name, hopts, copts, FLAGS_outdir, "dvz_central", "", "dV_{z}", "fraction");
//
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
