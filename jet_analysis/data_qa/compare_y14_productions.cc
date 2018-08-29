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
    TTreeReaderValue<int> ntrackswhft(reader, "ntrackswhft");
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
      evt.primwhft = *ntrackswhft;
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

std::set<std::pair<int, int>> GetMatchedKeys(std::vector<std::unique_ptr<event_map>>& maps) {
  std::set<std::pair<int, int>> ret;

  for (auto& event : *maps.front()) {

    auto& key = event.first;
    
    bool use_event = true;
    for (int i = 1; i < maps.size(); ++i) {
      if ((*maps[i]).find(key) == (*maps[i]).end()) {
        use_event = false;
        break;
      }
    }
    if (use_event) ret.insert(key);
  }
  return ret;
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
  
  auto input_files = ParseArgStringToVec<string>(FLAGS_data);
  auto identifiers = ParseArgStringToVec<string>(FLAGS_dataNames);
  
  if (input_files.size() == 0) {
    LOG(ERROR) << "No input files specified";
    return 1;
  }
  
  // check to make sure the input files exist
  for (auto& file : input_files) {
    if (!boost::filesystem::exists(file)) {
      LOG(ERROR) << "Input file does not exist: " << file;
      return 1;
    }
  }
  
  if (identifiers.size() != input_files.size()) {
    LOG(ERROR) << "Need a single identifier for each input file: "
               << identifiers.size() << " received, but "
               << input_files.size() << "input files exist";
    return 1;
  }
  
  
  // build output directory if it doesn't exist, using boost::filesystem
  if (FLAGS_outdir.empty())
    FLAGS_outdir = "tmp";
  boost::filesystem::path dir(FLAGS_outdir.c_str());
  boost::filesystem::create_directories(dir);
  
  // load in the datasets
  std::vector<std::unique_ptr<event_map>> data;
  data.reserve(input_files.size());
  for (auto& file : input_files) {
    LOG(INFO) << "Loading file: " << file;
    data.push_back(std::make_unique<event_map>());
    if (!ReadInTree(file, FLAGS_tree, *data.back())) {
      LOG(ERROR) << "Failed to read in " << file << ": exiting";
      return 1;
    }
    LOG(INFO) << "Loaded: " << data.back()->size() << " events";
  }
  LOG(INFO) << "Loading complete";
  
  // create centrality definition
  CentralityRun14 centrality;

  std::vector<TH1D*> h_vertices;
  std::vector<TH1D*> h_vpdvz;
  std::vector<TH1D*> h_zdcrate;
  std::vector<TH1D*> h_nglobal;
  std::vector<TH1D*> h_nprim;
  std::vector<TH1D*> h_nprim_ana;
  std::vector<TH1D*> h_vz;
  std::vector<TH1D*> h_refmult;
  std::vector<TH1D*> h_rank;
  std::vector<TH1D*> h_dvz;
  std::vector<TH1D*> h_pt;
  std::vector<TH1D*> h_dca;
  std::vector<TH1D*> h_nhit;
  std::vector<TH1D*> h_nhitpos;
  std::vector<TH1D*> h_phi;
  std::vector<TH1D*> h_eta;
  std::vector<TH2D*> h_pt_cent;
  std::vector<TH1D*> h_cent;
  std::vector<TH1D*> h_zerocent_nprim;


  for (auto& prefix : identifiers) {
    h_vertices.push_back(new TH1D(MakeString(prefix, "vertices").c_str(),
                                                ";N_{vertices}", 50, 0, 50));
    h_vpdvz.push_back(new TH1D(MakeString(prefix, "vpdvz").c_str(),
                                             ";VPD V_{z} [cm]", 70, -35, 35));
    h_zdcrate.push_back(new TH1D(MakeString(prefix, "zdcrate").c_str(),
                                               ";zdcX [kHz]", 100, 0, 50));
    h_nglobal.push_back(new TH1D(MakeString(prefix, "global").c_str(),
                                                ";N_{global}", 200, 0, 2000));
    h_nprim.push_back(new TH1D(MakeString(prefix, "nprim").c_str(),
                                             ";N_{prim}", 200, 0, 2000));
    h_nprim_ana.push_back(new TH1D(MakeString(prefix, "nprimana").c_str(),
                                   ";N_{prim}", 200, 0, 2000));
    h_vz.push_back(new TH1D(MakeString(prefix, "vz").c_str(),
                                          ";V_{z} [cm]", 70, -35, 35));
    h_refmult.push_back(new TH1D(MakeString(prefix, "refmult").c_str(),
                                               ";refmult", 200, 0, 800));
    h_rank.push_back(new TH1D(MakeString(prefix, "rank").c_str(),
                                            ";rank", 100, -10, 10));
    h_dvz.push_back(new TH1D(MakeString(prefix, "dvz").c_str(),
                                           ";dV_{z} [cm]", 100, -5, 5));
    h_pt.push_back(new TH1D(MakeString(prefix, "pt").c_str(),
                                          ";p_{T} [GeV/c]", 100, 0.0, 5.0));
    h_dca.push_back(new TH1D(MakeString(prefix, "dca").c_str(),
                                          ";DCA [cm]", 100, 0.0, 5.0));
    h_nhit.push_back(new TH1D(MakeString(prefix, "nhit").c_str(),
                                           ";nHit", 50, 0, 50));
    h_nhitpos.push_back(new TH1D(MakeString(prefix, "nhitpos").c_str(),
                                            ";nHitpos", 50, 0, 50));
    h_phi.push_back(new TH1D(MakeString(prefix, "phi").c_str(),
                             ";#phi", 50, TMath::Pi(), TMath::Pi()));
    h_eta.push_back(new TH1D(MakeString(prefix, "eta").c_str(),
                                           ";#eta", 50, -1.0, 1.0));
    h_pt_cent.push_back(new TH2D(MakeString(prefix, "ptcent").c_str(),
                                 ";cent;p_{T}", 9, 0, 9, 50, 0, 5.0));
    h_cent.push_back(new TH1D(MakeString(prefix, "cent").c_str(), ";centrality",
                              9, 0, 9));
    h_zerocent_nprim.push_back(new TH1D(MakeString(prefix, "centzero").c_str(), ";centrality",
                              100, 0, 1000));
  }

  // first, get all keys that are present in all three datasets
  LOG(INFO) << "Finding matched events";
  std::set<std::pair<int, int>> matched_keys = GetMatchedKeys(data);
  LOG(INFO) << matched_keys.size() << " shared events between all datasets";

  // now we can loop over all datasets with all matched keys and plot
  // relevant data - first think we do is make sure no HFT info is in any
  // of them
  LOG(INFO) << "Starting event loop";
  for (auto& key : matched_keys) {
    
    // first make sure no hft info is in this event
    bool use = true;
    for (int i = 0; i < data.size(); ++i) {
    
      auto& map = *data[i];
      if (map[key].pxl + map[key].ist + map[key].ssd > 0) {
        use = false;
        continue;
      }
    }
    
    if (!use)
      continue;
    
    
    // fill histograms
    for (int i = 0; i < data.size(); ++i) {
      
      auto& map = *data[i];
      
      if (map[key].rank < 0)
        continue;
    
      centrality.setEvent(map[key].runid, map[key].refmult, map[key].zdcrate, map[key].vz);
      int cent = centrality.centrality9();
      h_cent[i]->Fill(cent);
        
      
      h_vertices[i]->Fill(map[key].nvertices);
      h_vpdvz[i]->Fill(map[key].vpdvz);
      h_zdcrate[i]->Fill(map[key].zdcrate / 1000.0);
      h_nglobal[i]->Fill(map[key].nglobal);
      h_nprim[i]->Fill(map[key].nprim);
      h_vz[i]->Fill(map[key].vz);
      h_refmult[i]->Fill(map[key].refmult);
      h_rank[i]->Fill(map[key].rank);
      h_dvz[i]->Fill(map[key].dvz);

      if (map[key].primwhft > 0)
        LOG(INFO) << "WTF MAN";
      
      int nprim_ana = 0;
      for (int j = 0; j < map[key].pt.size(); ++j) {

        if (fabs(map[key].eta[j]) > 1.0) continue;
        if (map[key].nhit[j] < 20) continue;
        if ((double)map[key].nhit[j] / map[key].nhitspos[j] < 0.52) continue;
        if ((double)map[key].dca[j] > 1.0) continue;
        nprim_ana++;
        h_pt[i]->Fill(map[key].pt[j]);
        h_dca[i]->Fill(map[key].dca[j]);
        h_nhit[i]->Fill(map[key].nhit[j]);
        h_nhitpos[i]->Fill(map[key].nhitspos[j]);
        h_eta[i]->Fill(map[key].eta[j]);
        h_phi[i]->Fill(map[key].phi[j]);
        h_pt_cent[i]->Fill(cent, map[key].pt[j]);
      }
      h_nprim_ana[i]->Fill(nprim_ana);
      if (cent == 0)
        h_zerocent_nprim[i]->Fill(nprim_ana);
    }
  }
  
  // normalize histograms
  for (int i = 0; i < h_vertices.size(); ++i) {
    h_vertices[i]->Scale(1.0 / h_vertices[i]->Integral());
    h_vpdvz[i]->Scale(1.0 / h_vpdvz[i]->Integral());
    h_zdcrate[i]->Scale(1.0 / h_zdcrate[i]->Integral());
    h_nglobal[i]->Scale(1.0 / h_nglobal[i]->Integral());
    h_nprim[i]->Scale(1.0 / h_nprim[i]->Integral());
    h_nprim_ana[i]->Scale(1.0 / h_nprim_ana[i]->Integral());
    h_vz[i]->Scale(1.0 / h_vz[i]->Integral());
    h_refmult[i]->Scale(1.0 / h_refmult[i]->Integral());
    h_rank[i]->Scale(1.0 / h_rank[i]->Integral());
    h_dvz[i]->Scale(1.0 / h_dvz[i]->Integral());
    h_pt[i]->Scale(1.0 / h_vz[i]->GetEntries());
    h_dca[i]->Scale(1.0 / h_dca[i]->Integral());
    h_nhit[i]->Scale(1.0 / h_nhit[i]->Integral());
    h_nhitpos[i]->Scale(1.0 / h_nhitpos[i]->Integral());
    h_eta[i]->Scale(1.0 / h_eta[i]->Integral());
    h_phi[i]->Scale(1.0 / h_phi[i]->Integral());
    h_zerocent_nprim[i]->Scale(1.0 / h_zerocent_nprim[i]->Integral());
  }
  
  // print out
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
  
  Overlay1D(h_vertices, identifiers, hopts, copts, FLAGS_outdir, "nvertices",
            "", "N_{vertices}", "fraction", "");
  Overlay1D(h_vpdvz, identifiers, hopts, copts, FLAGS_outdir, "vpdvz",
            "", "VPD V_{z} [cm]", "fraction", "");
  Overlay1D(h_zdcrate, identifiers, hopts, cOptsTopLeftLeg, FLAGS_outdir, "zdcrate",
            "", "zdcX [kHz]", "fraction", "");
  Overlay1D(h_nglobal, identifiers, hopts, coptslogy, FLAGS_outdir, "nglobal",
            "", "N_{global}", "fraction", "");
  Overlay1D(h_nprim, identifiers, hopts, coptslogy, FLAGS_outdir, "nprimary",
            "", "N_{primary}", "fraction", "");
  Overlay1D(h_nprim_ana, identifiers, hopts, coptslogy, FLAGS_outdir, "nprimaryana",
            "", "N_{primary}", "fraction", "");
  Overlay1D(h_vz, identifiers, hopts, cOptsBottomLeftLegLogy, FLAGS_outdir, "vz",
            "", "V_{z} [cm]", "fraction", "");
  Overlay1D(h_refmult, identifiers, hopts, coptslogy, FLAGS_outdir, "refmult",
            "", "refmult", "fraction", "");
  Overlay1D(h_rank, identifiers, hopts, coptslogy, FLAGS_outdir, "rank",
            "", "rank", "fraction", "");
  Overlay1D(h_dvz, identifiers, hopts, copts, FLAGS_outdir, "dvz",
            "", "dV_{z} [cm]", "fraction", "");
  Overlay1D(h_pt, identifiers, hopts, coptslogy, FLAGS_outdir, "pt",
            "", "p_{T}", "fraction", "");
  Overlay1D(h_dca, identifiers, hopts, coptslogy, FLAGS_outdir, "dca",
            "", "DCA [cm]", "fraction", "");
  Overlay1D(h_nhit, identifiers, hopts, cOptsTopLeftLeg, FLAGS_outdir, "nhit",
            "", "nHitsFit", "fraction", "");
  Overlay1D(h_nhitpos, identifiers, hopts, cOptsTopLeftLeg, FLAGS_outdir, "nhitpos",
            "", "nHitsPoss", "fraction", "");
  Overlay1D(h_eta, identifiers, hopts, cOptsBottomLeg, FLAGS_outdir, "eta",
            "", "#eta", "fraction", "");
  Overlay1D(h_phi, identifiers, hopts, cOptsBottomLeg, FLAGS_outdir, "phi",
            "", "#phi", "fraction", "");
  Overlay1D(h_zerocent_nprim, identifiers, hopts, coptslogy, FLAGS_outdir, "zerocent",
            "", "N_{primary}", "fraction", "");
  
  std::vector<std::vector<TH1D*>> cent_projections;
  for (int i = 1; i <= h_pt_cent[0]->GetXaxis()->GetNbins(); ++i) {
    cent_projections.push_back(std::vector<TH1D*>());
    for (int j = 0; j < h_pt_cent.size(); ++j) {
      string tmpname = identifiers[j] + "_" + std::to_string(i-1) + "ptproj";
      TH1D* tmp = (TH1D*) h_pt_cent[j]->ProjectionY(tmpname.c_str(), i, i);
      tmp->Scale(1.0 / h_cent[j]->GetBinContent(i));
      cent_projections[i-1].push_back(tmp);
    }
  }
  
  LOG(INFO) << "integral: " << cent_projections[0][0]->Integral(21, 50);
  LOG(INFO) << "integral: " << cent_projections[0][1]->Integral(21, 50);
  LOG(INFO) << "integral: " << cent_projections[0][2]->Integral(21, 50);
  
  for (int i = 0; i < cent_projections.size(); ++i) {
    Overlay1D(cent_projections[i], identifiers, hopts, coptslogy, FLAGS_outdir, "cent" + std::to_string(i),
              "", "p_{T}", "p_{T}/N_{events}", "");
  }
  
  TH1D* tmp1 = (TH1D*) cent_projections[0][0]->Clone("tmp1");
  TH1D* tmp2 = (TH1D*) cent_projections[0][1]->Clone("tmp2");
  tmp1->Divide(cent_projections[0][2]);
  tmp2->Divide(cent_projections[0][2]);
  TCanvas c;
  tmp1->GetYaxis()->SetRangeUser(0.8, 1.2);
  tmp2->GetYaxis()->SetRangeUser(0.8, 1.2);
  tmp1->Draw();
  c.SaveAs("tmp1.pdf");
  tmp2->Draw();
  c.SaveAs("tmp2.pdf");
  
  return 0;
}
