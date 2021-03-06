// differential_aj_auau.cxx

#include <iostream>
#include <string>
#include <fstream>
#include <unordered_map>
#include <random>
#include <exception>

#include "jet_analysis/util/arg_helper.hh"
#include "jet_analysis/util/trigger_lookup.hh"
#include "jet_analysis/util/reader_util.hh"
#include "jet_analysis/util/string_util.hh"
#include "jet_analysis/util/vector_conversion.hh"
#include "jet_analysis/efficiency/run14_eff.hh"
#include "jet_analysis/dijet_worker/dijet_worker.hh"
#include "jet_analysis/centrality/centrality_run14.hh"

#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TH2.h"
#include "TH1.h"
#include "TProfile2D.h"

#include "TStarJetPicoReader.h"
#include "TStarJetPicoEvent.h"
#include "TStarJetPicoEventHeader.h"
#include "TStarJetPicoEventCuts.h"
#include "TStarJetPicoTrackCuts.h"
#include "TStarJetPicoTowerCuts.h"
#include "TStarJetVectorContainer.h"
#include "TStarJetVector.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

using std::string;
struct Options {
  string name        = "job"; /* output file name */
  int id             = 0;     /* job id */
  string input       = "";    /* root file/root file list*/
  string reader      = "";    /* settings file for primary reader */
  string out_dir     = "";    /* directory to save output in */
  string tow_list    = "";    /* list of hot towers to remove */
  string run_list    = "";    /* list of runs to remove */
  string triggers    = "";    /* triggers to consider (see trigger_lookup.hh) */
  string const_eta   = "";    /* constituent eta cut */
  string lc_pt       = "";    /* leading hard constituent pt cut */
  string sc_pt       = "";    /* subleading hard constituent pt cut */
  string lcm_pt      = "";    /* leading matched constituent pt cut */
  string scm_pt      = "";    /* subleading matched constituent pt cut */
  string lead_r      = "";    /* lead jet radii */
  string sub_r       = "";    /* sublead jet radii */
  string lj_pt       = "";    /* leading hard jet pt cut */
  string sj_pt       = "";    /* subleading hard jet pt cut */
  double dca         = 1.0;   /* track dca cut */
};

int main(int argc, char* argv[]) {
  
  // parse command line options
  Options opts;
  for (int i = 1; i < argc; ++i) {
    if (ParseStrFlag(string(argv[i]), "--name", &opts.name) ||
        ParseIntFlag(string(argv[i]), "--id", &opts.id) ||
        ParseStrFlag(string(argv[i]), "--input", &opts.input) ||
        ParseStrFlag(string(argv[i]), "--readerSetting", &opts.reader) ||
        ParseStrFlag(string(argv[i]), "--outDir", &opts.out_dir) ||
        ParseStrFlag(string(argv[i]), "--towList", &opts.tow_list) ||
        ParseStrFlag(string(argv[i]), "--runList", &opts.run_list) ||
        ParseStrFlag(string(argv[i]), "--triggers", &opts.triggers) ||
        ParseStrFlag(string(argv[i]), "--constEta", &opts.const_eta) ||
        ParseStrFlag(string(argv[i]), "--leadConstPt", &opts.lc_pt) ||
        ParseStrFlag(string(argv[i]), "--subConstPt", &opts.sc_pt) ||
        ParseStrFlag(string(argv[i]), "--leadConstPtMatch", &opts.lcm_pt) ||
        ParseStrFlag(string(argv[i]), "--subConstPtMatch", &opts.scm_pt) ||
        ParseStrFlag(string(argv[i]), "--leadR", &opts.lead_r) ||
        ParseStrFlag(string(argv[i]), "--subR", &opts.sub_r) ||
        ParseStrFlag(string(argv[i]), "--leadJetPt", &opts.lj_pt) ||
        ParseStrFlag(string(argv[i]), "--subJetPt", &opts.sj_pt) ||
        ParseFloatFlag(string(argv[i]), "--DCA", &opts.dca)) continue;
    std::cerr << "Unknown command line option: " << argv[i] << std::endl;
    return 1;
  }
  
  // reader & environment initialization
  // -----------------------------------
  
  // check to make sure the input file paths are sane
  if (!boost::filesystem::exists(opts.input)) {
    std::cerr << "input file does not exist: " << opts.input << std::endl;;
    return 1;
  }
  
  // first, build our input chain
  TChain* chain = NewChainFromInput(opts.input);
  
  // build output directory if it doesn't exist, using boost::filesystem
  if (opts.out_dir.empty())
    opts.out_dir = "tmp";
  boost::filesystem::path dir(opts.out_dir.c_str());
  boost::filesystem::create_directories(dir);
  
  // create output file from the given directory, name & id
  string outfile_name = opts.out_dir + "/" + opts.name + MakeString(opts.id) + ".root";
  TFile out(outfile_name.c_str(), "RECREATE");
  
  // initialize the reader
  TStarJetPicoReader* reader = new TStarJetPicoReader();
  if (!opts.reader.empty()) {
    InitReader(reader, chain, opts.reader, opts.tow_list, opts.run_list);
  }
  else {
    InitReaderWithDefaults(reader, chain, opts.tow_list, opts.run_list);
  }
  reader->GetTrackCuts()->SetDCACut(opts.dca);
  
  // get the trigger IDs that will be used
  std::set<unsigned> triggers = GetTriggerIDs(opts.triggers);
  
  std::cout << "taking triggers: " << opts.triggers << " for primary" << std::endl;
  std::cout << "trigger ids: ";
  for (auto i : triggers)
    std::cout << i << " ";
  std::cout << std::endl;
  
  // parse jetfinding variables
  // --------------------------
  
  // first, hard code the algorithm to be anti-kt
  std::set<fastjet::JetAlgorithm> alg{fastjet::antikt_algorithm};
  
  // constituent range
  std::set<double> const_eta              = ParseArgString<double>(opts.const_eta);
  
  // leading jet
  std::set<double> lead_const_hard_pt     = ParseArgString<double>(opts.lc_pt);
  std::set<double> lead_const_match_pt    = ParseArgString<double>(opts.lcm_pt);
  std::set<double> lead_R                 = ParseArgString<double>(opts.lead_r);
  std::set<double> lead_hard_pt           = ParseArgString<double>(opts.lj_pt);
  
  // subleading jet
  std::set<double> sublead_const_hard_pt  = ParseArgString<double>(opts.sc_pt);
  std::set<double> sublead_const_match_pt = ParseArgString<double>(opts.scm_pt);
  std::set<double> sublead_R              = ParseArgString<double>(opts.sub_r);
  std::set<double> sublead_hard_pt        = ParseArgString<double>(opts.sj_pt);
  
  // here we can initialize the worker
  std::cout << "initializing worker..." << std::endl;
  DijetWorker worker(alg, lead_hard_pt, lead_R, sublead_hard_pt, sublead_R,
                     lead_const_hard_pt, lead_const_match_pt, sublead_const_hard_pt,
                     sublead_const_match_pt, const_eta);
  worker.Initialize();
  
  std::cout << "worker initialized - number of dijet definitions: "
            << worker.Size() <<std::endl;
  
  std::set<std::string> keys = worker.Keys();
  for (auto key : keys)
    std::cout << key << std::endl;
  
  // create an output tree for each definition
  // -----------------------------------------
  
  std::unordered_map<std::string, std::shared_ptr<TTree>> trees;
  
  // and the necessary branches
  std::unordered_map<std::string, int> run_id_dict;
  std::unordered_map<std::string, int> event_id_dict;
  std::unordered_map<std::string, double> vz_dict;
  std::unordered_map<std::string, int> refmult_dict;
  std::unordered_map<std::string, int> grefmult_dict;
  std::unordered_map<std::string, double> refmultcorr_dict;
  std::unordered_map<std::string, double> grefmultcorr_dict;
  std::unordered_map<std::string, int> cent_dict;
  std::unordered_map<std::string, double> zdcrate_dict;
  std::unordered_map<std::string, double> reactionplane_dict;
  std::unordered_map<std::string, int> nglobal_dict;
  std::unordered_map<std::string, int> npart_dict;
  std::unordered_map<std::string, TLorentzVector> lead_hard_jet_dict;
  std::unordered_map<std::string, int> lead_hard_jet_nconst_dict;
  std::unordered_map<std::string, double> lead_hard_rho_dict;
  std::unordered_map<std::string, double> lead_hard_sigma_dict;
  std::unordered_map<std::string, double> lead_hard_area_dict;
  std::unordered_map<std::string, TLorentzVector> lead_match_jet_dict;
  std::unordered_map<std::string, int> lead_match_jet_nconst_dict;
  std::unordered_map<std::string, double> lead_match_rho_dict;
  std::unordered_map<std::string, double> lead_match_sigma_dict;
  std::unordered_map<std::string, double> lead_match_area_dict;
  std::unordered_map<std::string, TLorentzVector> sublead_hard_jet_dict;
  std::unordered_map<std::string, int> sublead_hard_jet_nconst_dict;
  std::unordered_map<std::string, double> sublead_hard_rho_dict;
  std::unordered_map<std::string, double> sublead_hard_sigma_dict;
  std::unordered_map<std::string, double> sublead_hard_area_dict;
  std::unordered_map<std::string, TLorentzVector> sublead_match_jet_dict;
  std::unordered_map<std::string, int> sublead_match_jet_nconst_dict;
  std::unordered_map<std::string, double> sublead_match_rho_dict;
  std::unordered_map<std::string, double> sublead_match_sigma_dict;
  std::unordered_map<std::string, double> sublead_match_area_dict;
  
  // fill the maps first, so that they don't decide to resize/move themselves
  // after branch creation...
  for (auto key : keys) {
    run_id_dict.insert({key, 0});
    event_id_dict.insert({key, 0});
    vz_dict.insert({key, 0});
    refmult_dict.insert({key, 0});
    grefmult_dict.insert({key, 0});
    refmultcorr_dict.insert({key, 0});
    grefmultcorr_dict.insert({key, 0});
    cent_dict.insert({key,0});
    zdcrate_dict.insert({key, 0});
    reactionplane_dict.insert({key, 0});
    nglobal_dict.insert({key, 0});
    npart_dict.insert({key, 0});
    lead_hard_jet_dict.insert({key, TLorentzVector()});
    lead_hard_jet_nconst_dict.insert({key, 0});
    lead_hard_rho_dict.insert({key, 0});
    lead_hard_sigma_dict.insert({key, 0});
    lead_hard_area_dict.insert({key, 0});
    lead_match_jet_dict.insert({key, TLorentzVector()});
    lead_match_jet_nconst_dict.insert({key, 0});
    lead_match_rho_dict.insert({key, 0});
    lead_match_sigma_dict.insert({key, 0});
    lead_match_area_dict.insert({key, 0});
    sublead_hard_jet_dict.insert({key, TLorentzVector()});
    sublead_hard_jet_nconst_dict.insert({key, 0});
    sublead_hard_rho_dict.insert({key, 0});
    sublead_hard_sigma_dict.insert({key, 0});
    sublead_hard_area_dict.insert({key, 0});
    sublead_match_jet_dict.insert({key, TLorentzVector()});
    sublead_match_jet_nconst_dict.insert({key, 0});
    sublead_match_rho_dict.insert({key, 0});
    sublead_match_sigma_dict.insert({key, 0});
    sublead_match_area_dict.insert({key, 0});
  }
  
  for (auto key : keys) {
    std::shared_ptr<TTree> tmp = std::make_shared<TTree>(key.c_str(), key.c_str());
    // create branches for the tree
    tmp->Branch("runid", &run_id_dict[key]);
    tmp->Branch("eventid", &event_id_dict[key]);
    tmp->Branch("vz", &vz_dict[key]);
    tmp->Branch("refmult", &refmult_dict[key]);
    tmp->Branch("grefmult", &grefmult_dict[key]);
    tmp->Branch("refmultcorr", &refmultcorr_dict[key]);
    tmp->Branch("grefmultcorr", &grefmultcorr_dict[key]);
    tmp->Branch("cent", &cent_dict[key]);
    tmp->Branch("zdcrate", &zdcrate_dict[key]);
    tmp->Branch("rp", &reactionplane_dict[key]);
    tmp->Branch("nglobal", &nglobal_dict[key]);
    tmp->Branch("npart", &npart_dict[key]);
    tmp->Branch("jl", &lead_hard_jet_dict[key]);
    tmp->Branch("js", &sublead_hard_jet_dict[key]);
    tmp->Branch("jlm", &lead_match_jet_dict[key]);
    tmp->Branch("jsm", &sublead_match_jet_dict[key]);
    tmp->Branch("jlconst", &lead_hard_jet_nconst_dict[key]);
    tmp->Branch("jlrho", &lead_hard_rho_dict[key]);
    tmp->Branch("jlsig", &lead_hard_sigma_dict[key]);
    tmp->Branch("jlarea", &lead_hard_area_dict[key]);
    tmp->Branch("jlmconst", &lead_match_jet_nconst_dict[key]);
    tmp->Branch("jlmrho", &lead_match_rho_dict[key]);
    tmp->Branch("jlmsig", &lead_match_sigma_dict[key]);
    tmp->Branch("jlmarea", &lead_match_area_dict[key]);
    tmp->Branch("jsconst", &sublead_hard_jet_nconst_dict[key]);
    tmp->Branch("jsrho", &sublead_hard_rho_dict[key]);
    tmp->Branch("jssig", &sublead_hard_sigma_dict[key]);
    tmp->Branch("jsarea", &sublead_hard_area_dict[key]);
    tmp->Branch("jsmconst", &sublead_match_jet_nconst_dict[key]);
    tmp->Branch("jsmrho", &sublead_match_rho_dict[key]);
    tmp->Branch("jsmsig", &sublead_match_sigma_dict[key]);
    tmp->Branch("jsmarea", &sublead_match_area_dict[key]);
    
    trees.insert({key, tmp});
  }
  
  // histograms
  // ----------
  
  std::unordered_map<string, TH1D*> lead_jet_count_dict;
  std::unordered_map<string, TH1D*> sublead_jet_count_dict;
  
  for (auto key : keys) {
    // create a unique histogram name for each key
    string lead_name = key + "_lead_count";
    string sublead_name = key + "_sublead_count";
    TH1D* lead_tmp = new TH1D(lead_name.c_str(), "count lead jets", 800, 0.5, 800.5);
    TH1D* sublead_tmp = new TH1D(sublead_name.c_str(), "count sublead jets", 800, 0.5, 800.5);
    
    lead_jet_count_dict.insert({key, lead_tmp});
    sublead_jet_count_dict.insert({key, sublead_tmp});
  }
  
  TH1D* h_ntracks = new TH1D("ntracks", "", 100, 0, 100);
  fastjet::Selector track_pt_min_selector_tmp = fastjet::SelectorPtMin(2.0) && fastjet::SelectorAbsRapMax(1.0);
  
  // initialize centrality definition
  CentralityRun14 centrality;
  
  // define a selector to reject low momentum tracks
  fastjet::Selector track_pt_min_selector = fastjet::SelectorPtMin(0.2) && fastjet::SelectorAbsRapMax(1.0);
  
  // start the analysis loop
  // -----------------------
  std::set<std::pair<int, int>> suspect_events{{8139067, 2890}, {8133016, 3377}, {8159021, 24232}, {8174095, 4791}, {8141106, 85038}};
  try {
    while (reader->NextEvent()) {
      // Print out reader status every 10 seconds
      reader->PrintStatus(10);
      
      // headers for convenience
      TStarJetPicoEventHeader* header = reader->GetEvent()->GetHeader();
  
      // check if event fired a trigger we will use
      if (triggers.size() != 0) {
        bool use_event = false;
        for (auto trigger : triggers)
          if (header->HasTriggerId(trigger))
            use_event = true;
        if (!use_event) continue;
      }
      
      // get centrality for the event
      centrality.setEvent(header->GetRunId(), header->GetReferenceMultiplicity(),
                          header->GetZdcCoincidenceRate(), header->GetPrimaryVertexZ());
      double refmultcorr = centrality.refMultCorr();
      double centrality_bin = centrality.centrality16();
      
      // get the vector container
      TStarJetVectorContainer<TStarJetVector>* container = reader->GetOutputContainer();
      std::vector<fastjet::PseudoJet> primary_particles;
      ConvertTStarJetVector(container, primary_particles);
      
      // select tracks above the minimum pt threshold
      primary_particles = track_pt_min_selector(primary_particles);
      
      h_ntracks->Fill(track_pt_min_selector_tmp(primary_particles).size());
      
      // run the worker
      auto& worker_out = worker.Run(primary_particles);
      
      // process any found di-jet pairs
      for (auto& result : worker_out) {
        std::string key = result.first;
        ClusterOutput& out = result.second;
        if (out.found_lead)
          lead_jet_count_dict[key]->Fill(header->GetReferenceMultiplicity());
        if (out.found_sublead)
          sublead_jet_count_dict[key]->Fill(header->GetReferenceMultiplicity());
        
        // now fill dijet results
        if (out.found_match) {
          // fill all branches for that key
          run_id_dict[key] = header->GetRunId();
          event_id_dict[key] = header->GetEventId();
          vz_dict[key] = header->GetPrimaryVertexZ();
          refmult_dict[key] = header->GetReferenceMultiplicity();
          grefmult_dict[key] = header->GetGReferenceMultiplicity();
          refmultcorr_dict[key] = refmultcorr;
          grefmultcorr_dict[key] = header->GetCorrectedGReferenceMultiplicity();
          cent_dict[key] = centrality_bin;
          zdcrate_dict[key] = header->GetZdcCoincidenceRate();
          reactionplane_dict[key] = header->GetReactionPlaneAngle();
          nglobal_dict[key] = header->GetNGlobalTracks();
          npart_dict[key] = primary_particles.size();
          
          // set the four jets
          lead_hard_jet_dict[key] = TLorentzVector(out.lead_hard.px(),
                                                   out.lead_hard.py(),
                                                   out.lead_hard.pz(),
                                                   out.lead_hard.E());
          lead_hard_jet_nconst_dict[key] = out.lead_hard.constituents().size();
          lead_hard_rho_dict[key] = out.lead_hard_rho;
          lead_hard_sigma_dict[key] = out.lead_hard_sigma;
          lead_hard_area_dict[key] = out.lead_hard.area();
          lead_match_jet_dict[key] = TLorentzVector(out.lead_match.px(),
                                                    out.lead_match.py(),
                                                    out.lead_match.pz(),
                                                    out.lead_match.E());
          lead_match_jet_nconst_dict[key] = out.lead_match.constituents().size();
          lead_match_rho_dict[key] = out.lead_match_rho;
          lead_match_sigma_dict[key] = out.lead_match_sigma;
          lead_match_area_dict[key] = out.lead_match.area();
          sublead_hard_jet_dict[key] = TLorentzVector(out.sublead_hard.px(),
                                                      out.sublead_hard.py(),
                                                      out.sublead_hard.pz(),
                                                      out.sublead_hard.E());
          sublead_hard_jet_nconst_dict[key] = out.sublead_hard.constituents().size();
          sublead_hard_rho_dict[key] = out.sublead_hard_rho;
          sublead_hard_sigma_dict[key] = out.sublead_hard_sigma;
          sublead_hard_area_dict[key] = out.sublead_hard.area();
          sublead_match_jet_dict[key] = TLorentzVector(out.sublead_match.px(),
                                                       out.sublead_match.py(),
                                                       out.sublead_match.pz(),
                                                       out.sublead_match.E());
          sublead_match_jet_nconst_dict[key] = out.sublead_match.constituents().size();
          sublead_match_rho_dict[key] = out.sublead_match_rho;
          sublead_match_sigma_dict[key] = out.sublead_match_sigma;
          sublead_match_area_dict[key] = out.sublead_match.area();
          
          trees[key]->Fill();
        }
      }
    }
  } catch(std::exception& e) {
    std::cerr << "Caught: " << e.what() << " during analysis loop." << std::endl;
  }
  
  out.Write();
  out.Close();
  
  
  return 0;
}
