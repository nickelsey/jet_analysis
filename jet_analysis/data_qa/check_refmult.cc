// check_refmult.cc
// QA to see if the default refmult is ok
// or needs to be recalculated

#include "jet_analysis/util/arg_helper.hh"
#include "jet_analysis/util/trigger_lookup.hh"
#include "jet_analysis/util/reader_util.hh"

#include <string>

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

#include "TFile.h"
#include "TH1.h"

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
  string triggers    = "";       /* triggers to consider (see trigger_lookup.hh) */
  int begin_run      = 0;        /* reject runs before this */
  int end_run        = 99999999; /* reject runs after this */
};

int main(int argc, char* argv[]) {
  
  // parse command line options
  Options opts;
  for (int i = 1; i < argc; ++i) {
    if (ParseStrFlag(string(argv[i]), "--name", &opts.name) ||
        ParseStrFlag(string(argv[i]), "--id", &opts.id) ||
        ParseStrFlag(string(argv[i]), "--input", &opts.input) ||
        ParseStrFlag(string(argv[i]), "--outDir", &opts.out_dir) ||
        ParseStrFlag(string(argv[i]), "--towList", &opts.tow_list) ||
        ParseStrFlag(string(argv[i]), "--runList", &opts.run_list) ||
        ParseStrFlag(string(argv[i]), "--triggers", &opts.triggers) ||
        ParseIntFlag(string(argv[i]), "--beginRun", &opts.begin_run) ||
        ParseIntFlag(string(argv[i]), "--endRun", &opts.end_run) ||) continue;
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
  
  // build our input chain
  TChain* chain = NewChainFromInput(opts.input);
  
  // create output file from the given directory, name & id
  string outfile_name = opts.out_dir + "/" + opts.name + opts.id + ".root";
  TFile out(outfile_name.c_str(), "RECREATE");
  
  // initialize the reader(s)
  TStarJetPicoReader* reader = new TStarJetPicoReader();
  InitReaderWithDefaults(reader, chain, opts.tow_list, opts.run_list);
  
  // but we're going to change some settings for track cuts, so that we get all tracks
  reader->GetTrackCuts()->SetDCACut(3.0);                // distance of closest approach to primary vtx
  reader->GetTrackCuts()->SetMinNFitPointsCut(10);       // minimum fit points in track reco
  reader->GetTrackCuts()->SetFitOverMaxPointsCut(0.1);  // minimum ratio of fit points used over possible
  reader->GetTrackCuts()->SetMaxPtCut(1000);             // essentially infinity - cut in eventcuts
  
  // get the triggers IDs that will be used
  std::set<unsigned> triggers = GetTriggerIDs(opts.triggers);
  
  // create histograms
  TH1D* h_refmult = new TH1D("refmult", "refmult;refmult;counts", 800, 0.5, 800.5);
  TH1D* h_grefmult = new TH1D("grefmult", "grefmult;grefmult;counts", 800, 0.5, 800.5);
  TH1D* h_refmult_10 = new TH1D("refmult10", "refmult10;refmult;counts", 800, 0.5, 800.5);
  TH1D* h_refmult_11 = new TH1D("refmult11", "refmult11;refmult;counts", 800, 0.5, 800.5);
  TH1D* h_refmult_12 = new TH1D("refmult12", "refmult12;refmult;counts", 800, 0.5, 800.5);
  TH1D* h_refmult_13 = new TH1D("refmult13", "refmult13;refmult;counts", 800, 0.5, 800.5);
  TH1D* h_refmult_14 = new TH1D("refmult14", "refmult14;refmult;counts", 800, 0.5, 800.5);
  TH1D* h_refmult_15 = new TH1D("refmult15", "refmult15;refmult;counts", 800, 0.5, 800.5);
  
  TH1D* h_refmult_vpd30 = new TH1D("refmultvpd30", "refmultvpd30;refmult;counts", 800, 0.5, 800.5);
  TH1D* h_refmult_vpd5 = new TH1D("refmultvpd5", "refmultvpd5;refmult;counts", 800, 0.5, 800.5);
  TH1D* h_refmult_mbmon = new TH1D("refmultmb", "refmultmb;refmult;counts", 800, 0.5, 800.5);
  
  
  try {
    while (reader->NextEvent()) {
      
      // Print out reader status every 10 seconds
      reader->PrintStatus(10);
      
      // headers for convenience
      TStarJetPicoEventHeader* header = reader->GetEvent()->GetHeader();
      
      if (header->GetRunNumber() < opts.begin_run ||
          header->GetRunNumber() > opts.end_run)
        continue;
      
      int refmult = header->GetReferenceMultiplicity();
      int grefmult = header->GetGReferenceMultiplicity();
      
      if (header->HasTriggerId(450008) ||
          header->HasTriggerId(450018))
        h_refmult_vpd5->Fill(refmult);
      
      if (header->HasTriggerId(450010) ||
          header->HasTriggerId(450020))
        h_refmult_vpd30->Fill(refmult);
      
      if (header->HasTriggerId(450011) ||
          header->HasTriggerId(450021))
        h_refmult_mbmon->Fill(refmult);
      
      // fill these two
      h_refmult->Fill(refmult);
      h_grefmult->Fill(grefmult);
      
      // create some counters
      int refmult_10 = 0;
      int refmult_11 = 0;
      int refmult_12 = 0;
      int refmult_13 = 0;
      int refmult_14 = 0;
      int refmult_15 = 0;
      
      // get the vector container
      // get tracks & towers
      TList* tracks = reader->GetListOfSelectedTracks();
      TIter nextTrack(tracks);
      while(TStarJetPicoPrimaryTrack* track = (TStarJetPicoPrimaryTrack*) nextTrack()){
        
        if(fabs(track->GetEta()) > 0.5)
          continue;
        int hits = track->GetNOfFittedHits();
        if (hits > 10)
          refmult_10++;
        if (hits > 11)
          refmult_11++;
        if (hits > 12)
          refmult_12++;
        if (hits > 13)
          refmult_13++;
        if (hits > 14)
          refmult_14++;
        if (hits > 15)
          refmult_15++;
      }
      
      h_refmult_10->Fill(refmult_10);
      h_refmult_11->Fill(refmult_11);
      h_refmult_12->Fill(refmult_12);
      h_refmult_13->Fill(refmult_13);
      h_refmult_14->Fill(refmult_14);
      h_refmult_15->Fill(refmult_15);
    }
  } catch(std::exception& e) {
    std::cerr << "Caught: " << e.what() << " during analysis loop." << std::endl;
  }
  
  out.Write();
  out.Close();
  
  return 0;
}
