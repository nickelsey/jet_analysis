// compare_auau.cc

#include "jet_analysis/util/common.hh"

#include "jet_analysis/util/trigger_lookup.hh"
#include "jet_analysis/util/reader_util.hh"
#include "jet_analysis/util/string_util.hh"
#include "jet_analysis/centrality/centrality_run14.hh"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"


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


DEFINE_bool(y14, false, "select year 14 centrality");
DEFINE_bool(y11, false, "select year 11 centrality");
DEFINE_bool(y7, false, "select year 7 centrality");

DEFINE_string(name, "job", "job name");
DEFINE_int32(id, 0, "job id");
DEFINE_string(input, "", "input root/txt file");
DEFINE_string(outdir, "tmp", "output directory");
DEFINE_string(towlist, "", "bad tower list");
DEFINE_string(runlist, "", "bad run list");
DEFINE_string(triggers, "", "trigger string (see util/trigger_lookup.hh");
DEFINE_int32(nhits, 10, "minimum number of fit points");
DEFINE_double(eta, 1.0, "maximum eta for tracks");
DEFINE_double(dca, 3.0, "max DCA for tracks");
DEFINE_double(fitfrac, 0.0, "minimum fit points/fit points possible");

// for year 7 centrality
int year7Centrality(int refmult) {
  static std::vector<std::pair<int, int>> CentBoundariesY7{ {485, 1000}, {399, 484}, {269, 398}, {178, 268}, {114, 177},
                                                           {69, 113}, {39, 68}, {21, 38}, {10, 20}};
  for (int i = 0; i < CentBoundariesY7.size(); ++i) {
    if (refmult >= CentBoundariesY7[i].first &&
        refmult < CentBoundariesY7[i].second)
      return i;
  }
  return -1;
}

int main(int argc, char* argv[]) {
  
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);
  
  if ((FLAGS_y7 == false && FLAGS_y11 == false && FLAGS_y14 == false) ||
      (FLAGS_y7 && FLAGS_y11) ||
      (FLAGS_y7 && FLAGS_y14) ||
      (FLAGS_y11 && FLAGS_y14)) {
    std::cerr << "need to select year 14, year 11, or year 7" << std::endl;
    return 1;
  }
  
  // check to make sure the input file exists
  if (!boost::filesystem::exists(FLAGS_input)) {
    std::cerr << "input file does not exist: " << FLAGS_input << std::endl;;
    return 1;
  }
  
  // build output directory if it doesn't exist, using boost::filesystem
  if (FLAGS_outdir.empty())
    FLAGS_outdir = "tmp";
  boost::filesystem::path dir(FLAGS_outdir.c_str());
  boost::filesystem::create_directories(dir);
  
  // build our input chain
  TChain* chain = NewChainFromInput(FLAGS_input);
  
  // create output file from the given directory, name & id
  string outfile_name = FLAGS_outdir + "/" + FLAGS_name + MakeString(FLAGS_id) + ".root";
  TFile out(outfile_name.c_str(), "RECREATE");
  
  // initialize the reader(s)
  TStarJetPicoReader* reader = new TStarJetPicoReader();
  InitReaderWithDefaults(reader, chain, FLAGS_towlist, FLAGS_runlist);
  reader->GetTrackCuts()->SetDCACut(FLAGS_dca);                   // distance of closest approach to primary vtx
  reader->GetTrackCuts()->SetMinNFitPointsCut(FLAGS_nhits);       // minimum fit points in track reco
  reader->GetTrackCuts()->SetFitOverMaxPointsCut(FLAGS_fitfrac);  // minimum ratio of fit points used over possible
  reader->GetTrackCuts()->SetMaxPtCut(1000);                      // essentially infinity - cut in eventcuts
  
  // get the triggers IDs that will be used
  std::set<unsigned> triggers = GetTriggerIDs(FLAGS_triggers);
  
  // centrality
  
  // initialize centrality definition for run 14
  // and set a lower limit on 0-20% for year 7
  CentralityRun14 centrality;
  
  // run 7 centrality is handled by the function Year7Centrality(refmult)
  // run 11 centrality is baked into the headers via refcent
  
  // Histograms will calculate gaussian errors
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();
  
  TH3D* fulltracks = new TH3D("tracks", ";nglobal;refmult;nprimary", 400, 0, 4000, 400, 0, 800, 400, 0, 1200);
  TH3D* lumiTracks = new TH3D("lumitracks", ";zdc rate [kHz];nglobal;nprimary", 100, 0, 100, 400, 0, 4000, 400, 0, 1200);
  TH2D* recalcRefMult = new TH2D("refmult", ";refmult;recalc refmult", 400, 0, 800, 400, 0, 800);
  
  TH3D* fulltracks_vz = new TH3D("tracksvz", ";nglobal;refmult;nprimary", 400, 0, 4000, 400, 0, 800, 400, 0, 1200);
  TH3D* lumiTracks_vz = new TH3D("lumitracksvz", ";zdc rate [kHz];nglobal;nprimary", 100, 0, 100, 400, 0, 4000, 400, 0, 1200);
  TH2D* recalcRefMult_vz = new TH2D("refmultvz", ";refmult;recalc refmult", 400, 0, 800, 400, 0, 800);
  
  TH1D* pt = new TH1D("pt", ";pt", 100, 0, 10);
  TH1D* nhit = new TH1D("nhit", ";nhit", 50, 0, 50);
  TH1D* nhitpos = new TH1D("nhitpos", ";nhitpos", 50, 0, 50);
  TH1D* fitfrac = new TH1D("fitfrac", ";fitfrac", 100, 0, 1.0);
  
  TH2D* etaphi = new TH2D("etaphi", ";eta;phi", 40, -1, 1, 40 -TMath::Pi(), TMath::Pi());
  
  
  // start the event loop
  // --------------------
  while(reader->NextEvent()) {
    
    // Print out reader status every 10 seconds
    reader->PrintStatus(10);
    
    TStarJetPicoEvent* event = reader->GetEvent();
    TStarJetPicoEventHeader* header = event->GetHeader();
    
    int cent_bin = 0;
    int refmult = 0;
    int nglobal = header->GetNGlobalTracks();
    double luminosity = header->GetZdcCoincidenceRate() / 1000.0;
    if (FLAGS_y7) {
      cent_bin = year7Centrality(header->GetGReferenceMultiplicity());
      refmult = header->GetGReferenceMultiplicity();
    }
    else if (FLAGS_y11) {
      cent_bin = 8 - header->GetReferenceCentrality();
      refmult = header->GetReferenceMultiplicity();
    }
    else if (FLAGS_y14) {
      centrality.setEvent(header->GetRunId(), header->GetReferenceMultiplicity(),
                          header->GetZdcCoincidenceRate(), header->GetPrimaryVertexZ());
      cent_bin = centrality.centrality9();
      refmult = header->GetReferenceMultiplicity();
    }
    
    if (cent_bin < 0 || cent_bin > 8)
      continue;
    
    // do track loop
    TList* tracks = reader->GetListOfSelectedTracks();
    TIter nextTrack(tracks);
    int recalc_refmult = 0;
    int nprimary = 0;
    while(TStarJetPicoPrimaryTrack* track = (TStarJetPicoPrimaryTrack*) nextTrack()) {
      if (fabs(track->GetEta()) < 0.5 && track->GetDCA() < 3.0 && track->GetNOfFittedHits() >= 10)
        recalc_refmult++;
      if (fabs(track->GetEta()) < 1.0 && track->GetDCA() < 3.0 && track->GetNOfFittedHits() >= 10)
        nprimary++;
      
      pt->Fill(track->GetPt());
      nhit->Fill(track->GetNOfFittedHits());
      nhitpos->Fill(track->GetNOfPossHits());
      fitfrac->Fill((double)track->GetNOfFittedHits() / track->GetNOfPossHits());
      etaphi->Fill(track->GetEta(), track->GetPhi());
      
    }
    recalcRefMult->Fill(refmult, recalc_refmult);
    fulltracks->Fill(nglobal, refmult, nprimary);
    lumiTracks->Fill(luminosity, refmult, nprimary);
    
    if (fabs(header->GetPrimaryVertexZ()) > 15.0) {
      recalcRefMult_vz->Fill(refmult, recalc_refmult);
      fulltracks_vz->Fill(nglobal, refmult, nprimary);
      lumiTracks_vz->Fill(luminosity, refmult, nprimary);
    }
    
  }
  
  out.Write();
  out.Close();
  
  gflags::ShutDownCommandLineFlags();
  return 0;
}
