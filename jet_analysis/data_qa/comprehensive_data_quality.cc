// comprehensive_data_quality.cc
// produces histograms pertaining to event, track
// and calorimeter health for a given run period

#include "jet_analysis/util/arg_helper.hh"
#include "jet_analysis/util/trigger_lookup.hh"
#include "jet_analysis/util/reader_util.hh"
#include "jet_analysis/util/string_util.hh"
#include "jet_analysis/efficiency/run14_eff.hh"
#include "jet_analysis/efficiency/run7_eff.hh"


#include <string>

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"

#include "TStarJetPicoReader.h"
#include "TStarJetPicoEvent.h"
#include "TStarJetPicoEventHeader.h"
#include "TStarJetPicoEventCuts.h"
#include "TStarJetPicoPrimaryTrack.h"
#include "TStarJetPicoTower.h"
#include "TStarJetPicoTrackCuts.h"
#include "TStarJetPicoTowerCuts.h"
#include "TStarJetVectorContainer.h"
#include "TStarJetVector.h"
#include "TStarJetPicoTriggerInfo.h"
#include "TStarJetPicoUtils.h"

using std::string;

struct Options {
  string name        = "job";    /* output file name */
  string id          = "0";      /* job id */
  string input       = "";       /* root file/root file list*/
  string out_dir     = "tmp";    /* directory to save output in */
  string tow_list    = "";       /* list of hot towers to remove */
  string run_list    = "";       /* list of runs to remove */
  string triggers    = "";       /* triggers to analyze (see trigger_lookup.hh) */
  string alt_trig    = "";       /* for more control, can pass a comma separated list of trigger IDs */
  string runid_in    = "";       /* root file containing ttree of runids */
  string lumi_str    = "";       /* prefix for histograms to separate different luminosities */
  string effcurves   = "";       /* file holding y14 efficiency curves */
  string y7effcurves = "";       /* file holding y7 efficiency curves */
  bool useY14Eff     = false;    /* use y14 efficiency curves when doing efficiency corrections */
  bool useY7Eff      = false;    /* use y7 efficiency curves when doing efficiency corrections */
};

// creates a sorted vector of unique runids from a
// tree, which is used to create 
std::vector<unsigned> PrepareRunIds(TFile* file) {
  
  std::vector<unsigned> runIds;
  
  TTree* tree = (TTree*) file->Get("runid");
  TBranch* branch = tree->GetBranch("runid");
  unsigned runId;
  branch->SetAddress(&runId);
  
  int entries = tree->GetEntries();
  for(int i = 0; i < entries; ++i) {
    tree->GetEntry(i);
    if(std::find(runIds.begin(), runIds.end(), runId) == runIds.end()) {
      runIds.push_back(runId);
    }
  }
  std::sort(runIds.begin(), runIds.end());
  
  return runIds;
}


int main(int argc, char* argv[]) {
  
  // parse command line options
  // --------------------------
  Options opts;
  for (int i = 1; i < argc; ++i) {
    if (ParseStrFlag(string(argv[i]), "--name", &opts.name) ||
        ParseStrFlag(string(argv[i]), "--id", &opts.id) ||
        ParseStrFlag(string(argv[i]), "--input", &opts.input) ||
        ParseStrFlag(string(argv[i]), "--outDir", &opts.out_dir) ||
        ParseStrFlag(string(argv[i]), "--towList", &opts.tow_list) ||
        ParseStrFlag(string(argv[i]), "--runList", &opts.run_list) ||
        ParseStrFlag(string(argv[i]), "--runIDs", &opts.runid_in) ||
        ParseStrFlag(string(argv[i]), "--triggers", &opts.triggers) ||
        ParseStrFlag(string(argv[i]), "--triggerIDs", &opts.alt_trig) ||
        ParseStrFlag(string(argv[i]), "--histPrefix", &opts.lumi_str) ||
        ParseStrFlag(string(argv[i]), "--y14Eff", &opts.effcurves) ||
        ParseStrFlag(string(argv[i]), "--y7Eff", &opts.y7effcurves)) continue;
    std::cerr << "Unknown command line option: " << argv[i] << std::endl;
    return 1;
  }
  
  if (opts.y7effcurves != "" && opts.effcurves != "") {
    std::cerr << "cant use both y7 and y14 efficiency curves: exiting" << std::endl;
    return 1;
  }
  if (opts.y7effcurves != "")
    opts.useY7Eff = true;
  else if (opts.effcurves != "")
    opts.useY14Eff = true;
  
  // initialization
  // --------------
  
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
  
  // build our input chain
  TChain* chain = NewChainFromInput(opts.input);
  
  // create output file from the given directory, name & id
  string outfile_name = opts.out_dir + "/" + opts.name + opts.id + ".root";
  TFile out(outfile_name.c_str(), "RECREATE");
  
  // initialize the reader(s)
  TStarJetPicoReader* reader = new TStarJetPicoReader();
  InitReaderWithDefaults(reader, chain, opts.tow_list, opts.run_list);
  
  // get the triggers IDs that will be used
  std::set<unsigned> triggers = GetTriggerIDs(opts.triggers);
  std::set<unsigned> alt_triggers = ParseArgString<unsigned>(opts.alt_trig);
  for (auto trig : alt_triggers)
    triggers.insert(trig);
  
  // sort runids
  TFile* runid_file = new TFile(opts.runid_in.c_str(), "READ");
  std::vector<unsigned> runids = PrepareRunIds(runid_file);
  int runID_bins         = runids.size();
  double runID_bin_width = 1.0;
  double runID_low_edge  = 0.5;
  double runID_high_edge = runID_low_edge + runID_bin_width * runID_bins;

  // load efficiency classes
  if (opts.useY7Eff && opts.useY14Eff) {
    std::cerr << "error: cant use both y14 and y7 efficiency corrections, exiting" << std::endl;
    return 1;
  }
  Run14Eff* run14Eff;
  if (opts.useY14Eff) {
    run14Eff = new Run14Eff();
    run14Eff->loadFile(opts.effcurves);
  }
  
  Run7Eff* run7Eff;
  if (opts.useY7Eff) {
    run7Eff = new Run7Eff();
    run7Eff->loadFile(opts.y7effcurves);
  }
  
  // define default y14 centrality cuts
  std::vector<int> y14_refcent_def{420, 364, 276, 212, 156, 108, 68, 44, 28, 12, 0};
  std::vector<int> y7_refcent_def{269};
  // define default y7 centrality cuts
  
  
  // define histograms
  // -----------------
  
  // change to output file
  out.cd();
  
  // Histograms will calculate gaussian errors
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();
  
  // histogram prefix to differentiate different data sets
  string prefix = opts.lumi_str;
  
  // reference multiplicity
  TH2D* runID_refmult = new TH2D(MakeString(prefix, "runidrefmult").c_str(), ";runID;refMult",
                                 runID_bins, runID_low_edge, runID_high_edge,
                                 800, 0, 800);
  TH2D* runID_grefmult = new TH2D(MakeString(prefix, "runidgrefmult").c_str(), ";runID;gRefMult",
                                  runID_bins, runID_low_edge, runID_high_edge,
                                  800, 0, 800);
  TH2D* runID_nprim = new TH2D(MakeString(prefix, "runidnprim").c_str(), ";runID;N_{primary}",
                               runID_bins, runID_low_edge, runID_high_edge,
                               800, 0, 1600);
  TH2D* zdc_refmult = new TH2D(MakeString(prefix, "zdcrefmult").c_str(), ";zdc[khz];refMult",
                               100, 0, 100,
                               800, 0, 800);
  TH2D* bbc_refmult = new TH2D(MakeString(prefix, "bbcrefmult").c_str(), ";BBC [kHz];refMult",
                               100, 0, 1000,
                               800, 0, 800);
  TH2D* zdc_grefmult = new TH2D(MakeString(prefix, "zdcgrefmult").c_str(), ";zdc[khz];gRefMult",
                                100, 0, 100,
                                800, 0, 800);
  TH2D* vz_refmult = new TH2D(MakeString(prefix, "vzrefmult").c_str(), ";V_{z};refMult",
                              100, -30, 30,
                              800, 0, 800);
  TH2D* ref_gref = new TH2D(MakeString(prefix, "refgrefmult").c_str(), ";gRefMult;refMult",
                            800, 0, 800,
                            800, 0, 800);
  TH2D* prim_glob = new TH2D(MakeString(prefix, "primglob").c_str(), ";N_{global};N_{primary}",
                             400, 0, 4000,
                             400, 0, 2000);
  
  // vertex
  TH2D* runID_vx = new TH2D(MakeString(prefix, "runidvx").c_str(), ";runID;V_{x}",
                            runID_bins, runID_low_edge, runID_high_edge,
                            100, -3, 3);
  TH2D* runID_vy = new TH2D(MakeString(prefix, "runidvy").c_str(), ";runID;V_{y}",
                            runID_bins, runID_low_edge, runID_high_edge,
                            100, -3, 3);
  TH2D* runID_vz = new TH2D(MakeString(prefix, "runidvz").c_str(), ";runID;V_{z}",
                            runID_bins, runID_low_edge, runID_high_edge,
                            140, -35, 35);
  TH2D* vx_vy = new TH2D(MakeString(prefix, "vxvy").c_str(), ";V_{x};V_{y}",
                         100, -3, 3,
                         100, -3, 3);
  TH2D* vx_vz = new TH2D(MakeString(prefix, "vxvz").c_str(), ";V_{x};V_{z}",
                         100, -3, 3,
                         140, -35, 35);
  TH2D* vy_vz = new TH2D(MakeString(prefix, "vyvz").c_str(), ";V_{y};V_{z}",
                         100, -3, 3,
                         140, -35, 35);
  TH2D* zdc_vz = new TH2D(MakeString(prefix, "zdcvz").c_str(), ";ZDC Rate[kHz];V_{z}[cm]",
                          100, 0, 100,
                          140, -35, 35);
  TH2D* zdc_vzvpdvz = new TH2D(MakeString(prefix, "vzvpdvz").c_str(), ";ZDC Rate[kHz];V_{z} - VPD V_{z}[cm]",
                               100, 0, 100,
                               100, -3, 3);
  
  // coincidence rate
  TH2D* runID_zdc = new TH2D(MakeString(prefix, "runidzdc").c_str(), ";runID;ZDC Rate[kHz]",
                             runID_bins, runID_low_edge, runID_high_edge,
                             100, 0, 100);
  TH2D* runID_bbc = new TH2D(MakeString(prefix, "runidbbc").c_str(), ";runID;BBC Rate[kHz]",
                             runID_bins, runID_low_edge, runID_high_edge,
                             100, 0, 1000);
  
  // tracks
  TH2D* runID_px = new TH2D(MakeString(prefix, "runidpx").c_str(), ";runID;p_{x}",
                            runID_bins, runID_low_edge, runID_high_edge,
                            200, -15, 15);
  TH2D* runID_py = new TH2D(MakeString(prefix, "runidpy").c_str(), ";runID;p_{y}",
                            runID_bins, runID_low_edge, runID_high_edge,
                            200, -15, 15);
  TH2D* runID_pz = new TH2D(MakeString(prefix, "runidpz").c_str(), ";runID;p_{z}",
                            runID_bins, runID_low_edge, runID_high_edge,
                            200, -15, 15);
  TH2D* runID_pt = new TH2D(MakeString(prefix, "runidpt").c_str(), ";runID;p_{t}",
                            runID_bins, runID_low_edge, runID_high_edge,
                            200, 0, 30);
  TH2D* px_py = new TH2D(MakeString(prefix, "pxpy").c_str(), ";p_{x};p_{y}",
                         200, -15, 15,
                         200, -15, 15);
  TH2D* zdc_px = new TH2D(MakeString(prefix, "zdcdpx").c_str(), ";ZDC Rate[kHz];p_{x}",
                          100, 0, 100,
                          200, -15, 15);
  TH2D* zdc_py = new TH2D(MakeString(prefix, "zdcpy").c_str(), ";ZDC Rate[kHz];p_{y}",
                          100, 0, 100,
                          200, -15, 15);
  TH2D* zdc_pz = new TH2D(MakeString(prefix, "zdcpz").c_str(), ";ZDC Rate[kHz];p_{z}",
                          100, 0, 100,
                          200, -15, 15);
  TH2D* zdc_pt = new TH2D(MakeString(prefix, "zdcpt").c_str(), ";ZDC Rate[kHz];p_{t}",
                          100, 0, 100,
                          200, 0, 30);
  TH2D* runID_dca = new TH2D(MakeString(prefix, "runiddca").c_str(), ";runID;DCA[cm]",
                             runID_bins, runID_low_edge, runID_high_edge,
                             50, 0, 3);
  TH2D* runID_fit = new TH2D(MakeString(prefix, "runidfit").c_str(), ";runID;fit points",
                             runID_bins, runID_low_edge, runID_high_edge,
                             50, 0, 50);
  TH2D* runID_tracketa = new TH2D(MakeString(prefix, "tracketa").c_str(), ";runID;#eta",
                                  runID_bins, runID_low_edge, runID_high_edge,
                                  100, -1, 1);
  TH2D* runID_trackphi = new TH2D(MakeString(prefix, "trackphi").c_str(), ";runID;#phi",
                                  runID_bins, runID_low_edge, runID_high_edge,
                                  100, -TMath::Pi(), TMath::Pi());
  TH2D* tracketa_phi = new TH2D(MakeString(prefix, "tracketaphi").c_str(), ";#eta;#phi",
                                100, -1.0, 1.0,
                                100, -TMath::Pi(), TMath::Pi());
  TH2D* ref_ptcorr = new TH2D(MakeString(prefix, "refmultptcorr").c_str(), ";refmult;p_{t}",
                              100, 0, 100,
                              200, 0, 30);
  
  // calorimeter
  TH2D* runID_towe = new TH2D(MakeString(prefix, "runidtowe").c_str(), ";runID;E",
                              runID_bins, runID_low_edge, runID_high_edge,
                              100, 0, 50);
  TH2D* runID_towet = new TH2D(MakeString(prefix, "runidtowet").c_str(), ";runID;E_{T}",
                               runID_bins, runID_low_edge, runID_high_edge,
                               100, 0, 50);
  TH2D* runID_towadc = new TH2D(MakeString(prefix, "runidtowadc").c_str(), ";runID;ADC",
                                runID_bins, runID_low_edge, runID_high_edge,
                                100, 0, 1000);
  TH2D* zdc_towe = new TH2D(MakeString(prefix, "zdctowe").c_str(), ";ZDC Rate[kHz];E",
                              100, 0, 100,
                              100, 0, 50);
  TH2D* zdc_towet = new TH2D(MakeString(prefix, "zdctowet").c_str(), ";ZDC Rate[kHz];E_{T}",
                               100, 0, 100,
                               100, 0, 50);
  TH2D* zdc_towadc = new TH2D(MakeString(prefix, "zdctowadc").c_str(), ";ZDC Rate[kHz];ADC",
                                100, 0, 100,
                                100, 0, 1000);
  TH2D* tow_towe = new TH2D(MakeString(prefix, "towtowe").c_str(), ";towerID;E",
                            4800, 0.5, 4800.5,
                            100, 0, 50);
  TH2D* tow_towet = new TH2D(MakeString(prefix, "towtowet").c_str(), ";towerID;E_{T}",
                             4800, 0.5, 4800.5,
                             100, 0, 50);
  TH2D* tow_towadc = new TH2D(MakeString(prefix, "towtowadc").c_str(), ";towerID;ADC",
                              4800, 0.5, 4800.5,
                              100, 0, 1000);
  TH2D* tow_etaphi = new TH2D(MakeString(prefix, "towetaphi").c_str(), ";#eta;#phi",
                              40, -1, 1,
                              120, -TMath::Pi(), TMath::Pi());
  
  // start the event loop
  // --------------------
  while(reader->NextEvent()) {
    // Print out reader status every 10 seconds
    reader->PrintStatus(10);
    
    // get the run ID & the map to the index we'll use for histograms
    unsigned runID = reader->GetEvent()->GetHeader()->GetRunId();
    int runidxmap = 0;
    if(std::find(runids.begin(), runids.end(), runID) != runids.end()) {
      runidxmap = std::distance(runids.begin(),  std::find(runids.begin(), runids.end(), runID));
    } else {
      std::cerr << "run id not listed, skipping event" << std::endl;
      continue;
    }
    TStarJetPicoEvent* event = reader->GetEvent();
    TStarJetPicoEventHeader* header = event->GetHeader();
   
    // check if event fired a trigger we will use
    if (triggers.size() != 0) {
      bool use_event = false;
      for (auto trigger : triggers)
        if (header->HasTriggerId(trigger))
          use_event = true;
      if (!use_event) continue;
    }
    
    // get tracks & towers
    TList* tracks = reader->GetListOfSelectedTracks();
    TIter nextTrack(tracks);
    TList* towers = reader->GetListOfSelectedTowers();
    TIter nextTower(towers);
    
    // get zdc rate in khz
    double zdc_khz = header->GetZdcCoincidenceRate()/1000.0;
    double bbc_khz = header->GetBbcCoincidenceRate()/1000.0;
    
    int refmult = header->GetReferenceMultiplicity();
    int cent = -1;
    if (opts.useY14Eff) {
      for (unsigned i = 0; i < y14_refcent_def.size(); ++i) {
        if (refmult >= y14_refcent_def[i]) {
          cent = i;
          break;
        }
      }
    }
    else if (opts.useY7Eff) {
      if (refmult > y7_refcent_def[0])
        cent = 0;
    }
    
    // fill event level histograms
    runID_refmult->Fill(runidxmap, header->GetReferenceMultiplicity());
    runID_grefmult->Fill(runidxmap, header->GetGReferenceMultiplicity());
    runID_nprim->Fill(runidxmap, header->GetNOfPrimaryTracks());
    zdc_refmult->Fill(zdc_khz, header->GetReferenceMultiplicity());
    bbc_refmult->Fill(bbc_khz, header->GetReferenceMultiplicity());
    zdc_grefmult->Fill(zdc_khz, header->GetGReferenceMultiplicity());
    vz_refmult->Fill(header->GetPrimaryVertexZ(), header->GetReferenceMultiplicity());
    ref_gref->Fill(header->GetReferenceMultiplicity(), header->GetGReferenceMultiplicity());
    prim_glob->Fill(header->GetNGlobalTracks(), header->GetNOfPrimaryTracks());
    runID_vx->Fill(runidxmap, header->GetPrimaryVertexX());
    runID_vy->Fill(runidxmap, header->GetPrimaryVertexY());
    runID_vz->Fill(runidxmap, header->GetPrimaryVertexZ());
    vx_vy->Fill(header->GetPrimaryVertexX(), header->GetPrimaryVertexY());
    vx_vz->Fill(header->GetPrimaryVertexX(), header->GetPrimaryVertexZ());
    vy_vz->Fill(header->GetPrimaryVertexY(), header->GetPrimaryVertexZ());
    zdc_vz->Fill(zdc_khz, header->GetPrimaryVertexZ());
    zdc_vzvpdvz->Fill(zdc_khz, header->GetPrimaryVertexZ() - header->GetVpdVz());
    runID_zdc->Fill(runidxmap, zdc_khz);
    runID_bbc->Fill(runidxmap, bbc_khz);
    
    while(TStarJetPicoPrimaryTrack* track = (TStarJetPicoPrimaryTrack*) nextTrack()) {
      runID_px->Fill(runidxmap, track->GetPx());
      runID_py->Fill(runidxmap, track->GetPy());
      runID_pz->Fill(runidxmap, track->GetPz());
      runID_pt->Fill(runidxmap, track->GetPt());
      px_py->Fill(track->GetPx(), track->GetPy());
      zdc_px->Fill(zdc_khz, track->GetPx());
      zdc_py->Fill(zdc_khz, track->GetPy());
      zdc_pz->Fill(zdc_khz, track->GetPz());
      zdc_pt->Fill(zdc_khz, track->GetPt());
      if (opts.useY14Eff) {
        ref_ptcorr->Fill(refmult, track->GetPt(), run14Eff->AuAuEff(track->GetPt(), track->GetEta(),
                                                                    cent,
                                                                    header->GetZdcCoincidenceRate()));
      }
      else if (opts.useY7Eff) {
        if (cent < 3)
        ref_ptcorr->Fill(refmult, track->GetPt(), run7Eff->AuAuEff020Avg(track->GetPt(), track->GetEta()));
      }
      runID_dca->Fill(runidxmap, track->GetDCA());
      runID_fit->Fill(runidxmap, track->GetNOfFittedHits());
      runID_tracketa->Fill(runidxmap, track->GetEta());
      runID_trackphi->Fill(runidxmap, track->GetPhi());
      tracketa_phi->Fill(track->GetEta(), track->GetPhi());
    }
    
    while(TStarJetPicoTower* tower = (TStarJetPicoTower*) nextTower()) {
      runID_towe->Fill(runidxmap, tower->GetEnergy());
      runID_towet->Fill(runidxmap, tower->GetEt());
      runID_towadc->Fill(runidxmap, tower->GetADC());
      zdc_towe->Fill(zdc_khz, tower->GetEnergy());
      zdc_towet->Fill(zdc_khz, tower->GetEt());
      zdc_towadc->Fill(zdc_khz, tower->GetADC());
      tow_towe->Fill(tower->GetId(), tower->GetEnergy());
      tow_towet->Fill(tower->GetId(), tower->GetEt());
      tow_towadc->Fill(tower->GetId(), tower->GetADC());
      tow_etaphi->Fill(tower->GetEta(), tower->GetPhi());
    }
    
  }
  
  out.Write();
  out.Close();
  return 0;
}

