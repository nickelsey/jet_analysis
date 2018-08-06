// p16id_p17id_compare.cc

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

#include "TCanvas.h"

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

// use glog and gflags
#include "gflags/gflags.h"
#include "glog/stl_logging.h"
#include "glog/logging.h"

// compare multiple productions of y14 on an
// event-by-event basis

DEFINE_string(name, "p16id_p17id", "output name");
DEFINE_string(data, "", "comma separated list of input root files");
DEFINE_string(dataNames, "", "a string identifier for each input file");
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
  double vz;
  int nprim;
  int refmult;
  double rank;
  double dvz;
  int primaries;
  int primwhft;
  std::vector<double> pt;
  std::vector<double> eta;
  std::vector<double> phi;
  std::vector<double> dca;
  std::vector<int> nhit;
  std::vector<int> nhitspos;
  std::vector<int> hft;
  
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
    TTreeReaderValue<double> vz(reader, "vz");
    TTreeReaderValue<int> refmult(reader, "refmult");
    TTreeReaderValue<int> nprim(reader, "nprim");
    TTreeReaderValue<double> rank(reader, "rank");
    TTreeReaderValue<double> dvz(reader, "dvz");
    TTreeReaderValue<double> ntrackswftf(reader, "ntrackswftf");
    TTreeReaderValue<std::vector<double>> dca(reader, "dca");
    TTreeReaderValue<std::vector<double>> pt(reader, "pt");
    TTreeReaderValue<std::vector<int>> nhits(reader, "nhits");
    TTreeReaderValue<std::vector<int>> nhitspos(reader, "nhitspos");
    TTreeReaderValue<std::vector<double>> eta(reader, "eta");
    TTreeReaderValue<std::vector<double>> phi(reader, "phi");
    TTreeReaderValue<std::vector<int>> hft(reader, "hft");
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
      evt.primwhft = *ntrackswftf;
      evt.dvz = *dvz;
      evt.dca = *dca;
      evt.pt = *pt;
      evt.nhit = *nhits;
      evt.nhitspos = *nhitspos;
      evt.eta = *eta;
      evt.phi = *phi;
      evt.hft = *hft;
      
      map[{*runid, *eventid}] = evt;
    }
    
    return true;
  }
  else
    return false;
}

int main(int argc, char* argv[]) {
  
  string usage = "compares productions of run 14 on an event-by-event basis";
  
  gflags::SetUsageMessage(usage);
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);
  
  // set drawing preferences for histograms and graphs
  gStyle->SetOptStat(false);
  gStyle->SetOptFit(false);
  gStyle->SetLegendBorderSize(0);
  
  auto input_files = ParseArgStringToVec<std::string>(FLAGS_data);
  auto identifiers = ParseArgStringToVec<std::string>(FLAGS_dataNames);
  
  if (input_files.size() == 0) {
    LOG(ERROR) << "no input files specified";
    return 1;
  }
  
  // check to make sure the input files exist
  for (auto& file : input_files) {
    if (!boost::filesystem::exists(file)) {
      LOG(ERROR) << "input file does not exist: " << file;
      return 1;
    }
  }
  
  if (identifiers.size() != input_files.size()) {
    LOG(ERROR) << "need a single identifier for each input file: "
               << identifiers.size() << " received, but "
               << input_files.size() << "input files exist";
    return 1;
  }
  
  
  // build output directory if it doesn't exist, using boost::filesystem
  if (FLAGS_outdir.empty())
    FLAGS_outdir = "tmp";
  boost::filesystem::path dir(FLAGS_outdir.c_str());
  boost::filesystem::create_directories(dir);
  
  return 0;
}
