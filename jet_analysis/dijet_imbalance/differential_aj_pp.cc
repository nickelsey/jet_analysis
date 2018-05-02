// differential_aj_pp.cxx

#include <iostream>
#include <string>
#include <fstream>
#include <unordered_map>
#include <random>
#include <exception>
#include <cmath>

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
#include "TStarJetPicoPrimaryTrack.h"
#include "TStarJetPicoTower.h"
#include "TStarJetPicoTrackCuts.h"
#include "TStarJetPicoTowerCuts.h"
#include "TStarJetVectorContainer.h"
#include "TStarJetVector.h"
#include "TStarJetPicoTriggerInfo.h"
#include "TStarJetPicoUtils.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/ClusterSequencePassiveArea.hh"
#include "fastjet/ClusterSequenceActiveArea.hh"
#include "fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/FunctionOfPseudoJet.hh"

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

using std::string;
struct Options {
  string name        = "job"; /* output file name */
  int id             = 0;     /* job id */
  string input       = "";    /* root file/root file list*/
  string embed       = "";    /* root file/file list for embedding*/
  string reader      = "";    /* settings file for primary reader */
  string readeremb   = "";    /* settings for embedding reader */
  string out_dir     = "";    /* directory to save output in */
  string tow_list    = "";    /* list of hot towers to remove */
  string run_list    = "";    /* list of runs to remove */
  string triggers    = "";    /* triggers to consider (see trigger_lookup.hh) */
  string trigemb     = "";    /* triggers to consider for embedding (should be mb) */
  string const_eta   = "";    /* constituent eta cut */
  string lc_pt       = "";    /* leading hard constituent pt cut */
  string sc_pt       = "";    /* subleading hard constituent pt cut */
  string lcm_pt      = "";    /* leading matched constituent pt cut */
  string scm_pt      = "";    /* subleading matched constituent pt cut */
  string lead_r      = "";    /* lead jet radii */
  string sub_r       = "";    /* sublead jet radii */
  string lj_pt       = "";    /* leading hard jet pt cut */
  string sj_pt       = "";    /* subleading hard jet pt cut */
  string effic_file  = "";    /* file containing efficiency curves for run 14 */
  bool use_effic     = true;  /* turn on efficiency corrections for pp data */
  int tower_unc      = 0;     /* flag for tower systematic uncertainty (off = 0, or 1, -1) */
  int track_unc      = 0;     /* flag for tracking systematic uncertainty (off = 0, or 1, -1) */
  int embed_tries    = 32;    /* number of times we attempt to embed a trigger event */
};

int main(int argc, char* argv[]) {
  
  // parse command line options
  Options opts;
  for (int i = 1; i < argc; ++i) {
    if (ParseStrFlag(string(argv[i]), "--name", &opts.name) ||
        ParseIntFlag(string(argv[i]), "--id", &opts.id) ||
        ParseStrFlag(string(argv[i]), "--input", &opts.input) ||
        ParseStrFlag(string(argv[i]), "--embed", &opts.embed) ||
        ParseStrFlag(string(argv[i]), "--readerSetting", &opts.reader) ||
        ParseStrFlag(string(argv[i]), "--embedReaderSetting", &opts.readeremb) ||
        ParseStrFlag(string(argv[i]), "--outDir", &opts.out_dir) ||
        ParseStrFlag(string(argv[i]), "--towList", &opts.tow_list) ||
        ParseStrFlag(string(argv[i]), "--runList", &opts.run_list) ||
        ParseStrFlag(string(argv[i]), "--triggers", &opts.triggers) ||
        ParseStrFlag(string(argv[i]), "--embedTriggers", &opts.trigemb) ||
        ParseStrFlag(string(argv[i]), "--efficFile", &opts.effic_file) ||
        ParseBoolFlag(string(argv[i]), "--efficiency", &opts.use_effic) ||
        ParseStrFlag(string(argv[i]), "--constEta", &opts.const_eta) ||
        ParseStrFlag(string(argv[i]), "--leadConstPt", &opts.lc_pt) ||
        ParseStrFlag(string(argv[i]), "--subConstPt", &opts.sc_pt) ||
        ParseStrFlag(string(argv[i]), "--leadConstPtMatch", &opts.lcm_pt) ||
        ParseStrFlag(string(argv[i]), "--subConstPtMatch", &opts.scm_pt) ||
        ParseStrFlag(string(argv[i]), "--leadR", &opts.lead_r) ||
        ParseStrFlag(string(argv[i]), "--subR", &opts.sub_r) ||
        ParseStrFlag(string(argv[i]), "--leadJetPt", &opts.lj_pt) ||
        ParseStrFlag(string(argv[i]), "--subJetPt", &opts.sj_pt) ||
        ParseIntFlag(string(argv[i]), "--towerUnc", &opts.tower_unc) ||
        ParseIntFlag(string(argv[i]), "--trackUnc", &opts.track_unc) ||
        ParseIntFlag(string(argv[i]), "--embedTries", &opts.embed_tries)) continue;
    std::cerr << "Unknown command line option: " << argv[i] << std::endl;
    return 1;
  }
  
  // one sanity check - our systematic uncertainties and efficiency ratios only
  // make sense if we're embedding data, so return if settings are mismatched
  if (opts.embed.empty() && (opts.use_effic || opts.track_unc || opts.tower_unc)) {
    std::cerr << "can't use efficiency ratios without embedding: exiting" << std::endl;
    return 1;
  }
  
  // if we are not embedding, we can't embed multiple times :)
  if (opts.embed.empty())
    opts.embed_tries = 1;
  
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
  
  // setup embedding if requested
  TChain* embed_chain = nullptr;
  TStarJetPicoReader* embed_reader = new TStarJetPicoReader();
  if (!opts.embed.empty()) {
    embed_chain = new TChain("JetTree");
    if (HasEnding(opts.embed, ".root")) {
      embed_chain->Add(opts.embed.c_str());
    }
    else if (HasEnding(opts.embed, ".list") || HasEnding(opts.embed, ".txt")) {
      // read in our possible embedding events
      std::vector<string> embed_files;
      std::ifstream embed_stream(opts.embed);
      string line;
      while (std::getline(embed_stream, line)) {
        embed_files.push_back(line);
      }
      
      // we want to have a good sample of events, about 40x more than the pp trigger
      int min_events = chain->GetEntries() * 40;
      
      // we need to randomly sample the list using the index as the seed so its reproducible
      std::mt19937 gen(opts.id);
      std::uniform_int_distribution<int> uniform_dist(0, embed_files.size());
      
      while (embed_files.size() && embed_chain->GetEntries() < min_events) {
        int idx = uniform_dist(gen);
        embed_chain->Add(embed_files[idx].c_str());
        embed_files.erase(embed_files.begin() + idx);
        uniform_dist = std::uniform_int_distribution<>(0, embed_files.size());
      }
    }
    else {
      std::cerr << "error: unrecognized extension for embedding file type, exiting" << std::endl;
      return 1;
    }
    
    if (!opts.readeremb.empty()) {
      InitReader(embed_reader, embed_chain, opts.readeremb, opts.tow_list, opts.run_list);
    }
    else {
      InitReaderWithDefaults(embed_reader, embed_chain, opts.tow_list, opts.run_list);
    }
  }
  
  // get the trigger IDs that will be used
  std::set<unsigned> triggers = GetTriggerIDs(opts.triggers);
  
  std::set<unsigned> embed_triggers;
  if (!opts.trigemb.empty())
    embed_triggers = GetTriggerIDs(opts.trigemb);
  
  std::cout << "taking triggers: " << opts.triggers << " for primary" << std::endl;
  std::cout << "trigger ids: ";
  for (auto i : triggers)
    std::cout << i << " ";
  std::cout << std::endl;
  if (embed_reader != nullptr) {
    std::cout << "taking triggers: " << opts.trigemb << " for embedding" << std::endl;
    std::cout << "trigger ids: ";
    for (auto i : embed_triggers)
      std::cout << i << " ";
    std::cout << std::endl;
  }
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
  std::unordered_map<std::string, TLorentzVector> lead_match_jet_dict;
  std::unordered_map<std::string, int> lead_match_jet_nconst_dict;
  std::unordered_map<std::string, double> lead_match_rho_dict;
  std::unordered_map<std::string, double> lead_match_sigma_dict;
  std::unordered_map<std::string, TLorentzVector> sublead_hard_jet_dict;
  std::unordered_map<std::string, int> sublead_hard_jet_nconst_dict;
  std::unordered_map<std::string, double> sublead_hard_rho_dict;
  std::unordered_map<std::string, double> sublead_hard_sigma_dict;
  std::unordered_map<std::string, TLorentzVector> sublead_match_jet_dict;
  std::unordered_map<std::string, int> sublead_match_jet_nconst_dict;
  std::unordered_map<std::string, double> sublead_match_rho_dict;
  std::unordered_map<std::string, double> sublead_match_sigma_dict;
  
  // if embedding is being done
  std::unordered_map<std::string, int> embed_runid_dict;
  std::unordered_map<std::string, int> embed_eventid_dict;
  std::unordered_map<std::string, int> embed_refmult_dict;
  std::unordered_map<std::string, int> embed_grefmult_dict;
  std::unordered_map<std::string, double> embed_refmultcorr_dict;
  std::unordered_map<std::string, double> embed_grefmultcorr_dict;
  std::unordered_map<std::string, int> embed_cent_dict;
  std::unordered_map<std::string, int> embed_npart_dict;
  std::unordered_map<std::string, double> embed_vz_dict;
  std::unordered_map<std::string, double> embed_zdc_dict;
  std::unordered_map<std::string, double> embed_rp_dict;
  
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
    lead_match_jet_dict.insert({key, TLorentzVector()});
    lead_match_jet_nconst_dict.insert({key, 0});
    lead_match_rho_dict.insert({key, 0});
    lead_match_sigma_dict.insert({key, 0});
    sublead_hard_jet_dict.insert({key, TLorentzVector()});
    sublead_hard_jet_nconst_dict.insert({key, 0});
    sublead_hard_rho_dict.insert({key, 0});
    sublead_hard_sigma_dict.insert({key, 0});
    sublead_match_jet_dict.insert({key, TLorentzVector()});
    sublead_match_jet_nconst_dict.insert({key, 0});
    sublead_match_rho_dict.insert({key, 0});
    sublead_match_sigma_dict.insert({key, 0});
    
    if (embed_reader != nullptr) {
      embed_runid_dict.insert({key, 0});
      embed_eventid_dict.insert({key, 0});
      embed_vz_dict.insert({key, 0});
      embed_zdc_dict.insert({key, 0});
      embed_rp_dict.insert({key, 0});
      embed_cent_dict.insert({key, 0});
      embed_refmult_dict.insert({key, 0});
      embed_grefmult_dict.insert({key, 0});
      embed_refmultcorr_dict.insert({key, 0});
      embed_grefmultcorr_dict.insert({key, 0});
      embed_npart_dict.insert({key, 0});
    }
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
    tmp->Branch("jlmconst", &lead_match_jet_nconst_dict[key]);
    tmp->Branch("jlmrho", &lead_match_rho_dict[key]);
    tmp->Branch("jlmsig", &lead_match_sigma_dict[key]);
    tmp->Branch("jsconst", &sublead_hard_jet_nconst_dict[key]);
    tmp->Branch("jsrho", &sublead_hard_rho_dict[key]);
    tmp->Branch("jssig", &sublead_hard_sigma_dict[key]);
    tmp->Branch("jsmconst", &sublead_match_jet_nconst_dict[key]);
    tmp->Branch("jsmrho", &sublead_match_rho_dict[key]);
    tmp->Branch("jsmsig", &sublead_match_sigma_dict[key]);
    
    if (embed_reader != nullptr) {
      tmp->Branch("embed_eventid", &embed_eventid_dict[key]);
      tmp->Branch("embed_runid", &embed_runid_dict[key]);
      tmp->Branch("embed_refmult", &embed_refmult_dict[key]);
      tmp->Branch("embed_grefmult", &embed_grefmult_dict[key]);
      tmp->Branch("embed_refmultcorr", &embed_refmultcorr_dict[key]);
      tmp->Branch("embed_grefmultcorr", &embed_grefmultcorr_dict[key]);
      tmp->Branch("embed_cent", &embed_cent_dict[key]);
      tmp->Branch("embed_npart", &embed_npart_dict[key]);
      tmp->Branch("embed_rp", &embed_rp_dict[key]);
      tmp->Branch("embed_zdcrate", &embed_zdc_dict[key]);
      tmp->Branch("embed_vz", &embed_vz_dict[key]);
    }
    
    trees.insert({key, tmp});
  }
  
  // histograms
  // ----------
  
  std::unordered_map<string, TH1D*> lead_jet_count_dict;
  std::unordered_map<string, TH1D*> sublead_jet_count_dict;
  
  TProfile2D* eff_ratio  = new TProfile2D("pp_eff_ratio",
                                          "average pp efficiency ratio;p_{T};#eta;ratio",
                                          200, 0, 5, 10, -1.0, 1.0);
  TH1D* frac_finite = new TH1D("finite", "", 20, 0, 1);
  
  for (auto key : keys) {
    // create a unique histogram name for each key
    string lead_name = key + "_lead_count";
    string sublead_name = key + "_sublead_count";
    TH1D* lead_tmp = new TH1D(lead_name.c_str(), "count lead jets", 800, 0.5, 800.5);
    TH1D* sublead_tmp = new TH1D(sublead_name.c_str(), "count sublead jets", 800, 0.5, 800.5);
    
    lead_jet_count_dict.insert({key, lead_tmp});
    sublead_jet_count_dict.insert({key, sublead_tmp});
  }
  
  // initialize centrality definition
  CentralityRun14 centrality;
  
  // initialize efficiency curves
  Run14Eff* efficiency;
  if (opts.effic_file.empty())
    efficiency = new Run14Eff();
  else
    efficiency = new Run14Eff(opts.effic_file);
  
  switch(opts.track_unc) {
    case 0 :
      efficiency->setSystematicUncertainty(TrackingUnc::NONE);
      break;
    case 1 :
      efficiency->setSystematicUncertainty(TrackingUnc::POSITIVE);
      break;
    case -1 :
      efficiency->setSystematicUncertainty(TrackingUnc::NEGATIVE);
      break;
    default:
      std::cerr << "undefined tracking efficiency setting, exiting" << std::endl;
      return 1;
  }
  
  // define our tower uncertainty scaling as well
  double tower_scale_percent = 0.02;
  if (opts.tower_unc != 1 && opts.tower_unc != 0 && opts.tower_unc != -1) {
    std::cerr << "undefined tower efficiency setting, exiting" << std::endl;
  }
  double tower_scale = 1.0 + tower_scale_percent * opts.tower_unc;
  
  // and we'll need a random number generator for randomly throwing away particles
  std::mt19937 gen(opts.id);
  std::uniform_real_distribution<> dis(0.0,1.0);
  
  // define a selector to reject low momentum tracks
  fastjet::Selector track_pt_min_selector = fastjet::SelectorPtMin(0.2) && fastjet::SelectorAbsRapMax(1.0);
  
  // start the analysis loop
  // -----------------------
  try {
    while (reader->NextEvent()) {
      // Print out reader status every 10 seconds
      reader->PrintStatus(10);
      
      // headers for convenience
      TStarJetPicoEventHeader* header = reader->GetEvent()->GetHeader();
      TStarJetPicoEventHeader* embed_header = nullptr;
      if (embed_reader != nullptr) {
        embed_header = embed_reader->GetEvent()->GetHeader();
      }
      
      // check if event fired a trigger we will use
      if (triggers.size() != 0) {
        bool use_event = false;
        for (auto trigger : triggers)
          if (header->HasTriggerId(trigger))
            use_event = true;
        if (!use_event) continue;
      }
      
      // so we will use this trigger event: we now loop over 2xcentrality bins (32) events, and use as
      // as many of them as we can, without ending up in the same centrality bin twice
      std::set<int> cent_bins_used;
      for (int i = 0; i < opts.embed_tries; ++i) {
        if (cent_bins_used.size() == 16)
          break;
        
        int refmult = header->GetReferenceMultiplicity();
        double refmultcorr = refmult;
        int centrality_bin = -1;
        int embed_centrality = -1;
        
        std::vector<fastjet::PseudoJet> particles;
        std::vector<fastjet::PseudoJet> embed_particles;
        std::vector<fastjet::PseudoJet> primary_particles;
        
        
        // if embedding, read next event, loop to event 0 if needed
        if (embed_reader) {
          if (!GetNextValidEvent(embed_reader, embed_triggers)) {
            std::cerr << "no events found for embedding, given trigger requirements: exiting" << std::endl;
            return 1;
          }
        
          // effective reference multiplicity is refmult of pp + refmult of embedding event
          refmult = header->GetReferenceMultiplicity() + embed_header->GetReferenceMultiplicity();
         
          // get centrality for the event, based on the embedded event's information and a
          // combined refmult
          centrality.setEvent(embed_header->GetRunId(), refmult,
                              embed_header->GetZdcCoincidenceRate(), embed_header->GetPrimaryVertexZ());
          refmultcorr = centrality.refMultCorr();
          centrality_bin = centrality.centrality16();
          
          // we'll also get the embedded event's centrality to be complete
          centrality.setEvent(embed_header->GetRunId(), embed_header->GetReferenceMultiplicity(),
                              embed_header->GetZdcCoincidenceRate(), embed_header->GetPrimaryVertexZ());
          embed_centrality = centrality.centrality16();
          
          // check to see if the centrality has been used before, or if its 80-100%
          if (centrality_bin < 0 || cent_bins_used.find(centrality_bin) != cent_bins_used.end())
            continue;
          
          // now get the embedding event first
          // get the vector container for embedding
          TStarJetVectorContainer<TStarJetVector>* container_embed = embed_reader->GetOutputContainer();
          ConvertTStarJetVector(container_embed, particles);
          ConvertTStarJetVector(container_embed, embed_particles);
        }
        
        // and now convert the pp - if there is any efficiency curves to apply, do it now
        TStarJetVectorContainer<TStarJetVector>* container = reader->GetOutputContainer();
        if (opts.use_effic) {
          double norm = 0;
          double counter = 0;
          TStarJetVector* sv;
          for (int i = 0; i < container->GetEntries(); ++i) {
            sv = container->Get(i);
           
            if (sv->GetCharge()) {
              double ratio = efficiency->ratio(sv->Pt(), sv->Eta(), centrality_bin,
                                               embed_header->GetZdcCoincidenceRate());
               norm++
              if (!std::isfinite(ratio))
                continue;
              counter++;
              eff_ratio->Fill(sv->Pt(), sv->Eta(), ratio);
              double random_ = dis(gen);
              
              if ( random_ > ratio ) {
                continue;
              }
            }
            fastjet::PseudoJet tmpPJ = fastjet::PseudoJet(*sv);
            if (sv->GetCharge() == 0)
              tmpPJ *= tower_scale;
            tmpPJ.set_user_index( sv->GetCharge() );
            particles.push_back(tmpPJ);
            primary_particles.push_back(tmpPJ);
          }
          finite->Fill(counter/norm);
        }
        // no efficiency, convert all
        else {
          ConvertTStarJetVector(container, primary_particles);
          particles.insert(particles.end(), primary_particles.begin(), primary_particles.end());
        }
        
        // select tracks above the minimum pt threshold
        particles = track_pt_min_selector(particles);
        
        // run the worker
        auto worker_out = worker.Run(particles);
        
        // process any found di-jet pairs
        for (auto result : worker_out) {
          std::string key = result.first;
          ClusterOutput out = result.second;
          
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
            
            // if embedding is being done
            if (embed_reader != nullptr) {
              embed_runid_dict[key] = embed_header->GetRunId();
              embed_eventid_dict[key] = embed_header->GetEventId();
              embed_refmult_dict[key] = embed_header->GetReferenceMultiplicity();
              embed_grefmult_dict[key] = embed_header->GetGReferenceMultiplicity();
              embed_refmultcorr_dict[key] = embed_header->GetCorrectedReferenceMultiplicity();
              embed_grefmultcorr_dict[key] = embed_header->GetCorrectedGReferenceMultiplicity();
              embed_cent_dict[key] = embed_centrality;
              embed_vz_dict[key] = embed_header->GetPrimaryVertexZ();
              embed_zdc_dict[key] = embed_header->GetZdcCoincidenceRate();
              embed_rp_dict[key] = embed_header->GetReactionPlaneAngle();
              embed_npart_dict[key] = embed_particles.size();
            }
            
            // set the four jets
            lead_hard_jet_dict[key] = TLorentzVector(out.lead_hard.px(),
                                                     out.lead_hard.py(),
                                                     out.lead_hard.pz(),
                                                     out.lead_hard.E());
            lead_hard_jet_nconst_dict[key] = out.lead_hard.constituents().size();
            lead_hard_rho_dict[key] = out.lead_hard_rho;
            lead_hard_sigma_dict[key] = out.lead_hard_sigma;
            lead_match_jet_dict[key] = TLorentzVector(out.lead_match.px(),
                                                      out.lead_match.py(),
                                                      out.lead_match.pz(),
                                                      out.lead_match.E());
            lead_match_jet_nconst_dict[key] = out.lead_match.constituents().size();
            lead_match_rho_dict[key] = out.lead_match_rho;
            lead_match_sigma_dict[key] = out.lead_match_sigma;
            sublead_hard_jet_dict[key] = TLorentzVector(out.sublead_hard.px(),
                                                        out.sublead_hard.py(),
                                                        out.sublead_hard.pz(),
                                                        out.sublead_hard.E());
            sublead_hard_jet_nconst_dict[key] = out.sublead_hard.constituents().size();
            sublead_hard_rho_dict[key] = out.sublead_hard_rho;
            sublead_hard_sigma_dict[key] = out.sublead_hard_sigma;
            sublead_match_jet_dict[key] = TLorentzVector(out.sublead_match.px(),
                                                         out.sublead_match.py(),
                                                         out.sublead_match.pz(),
                                                         out.sublead_match.E());
            sublead_match_jet_nconst_dict[key] = out.sublead_match.constituents().size();
            sublead_match_rho_dict[key] = out.sublead_match_rho;
            sublead_match_sigma_dict[key] = out.sublead_match_sigma;
            
            trees[key]->Fill();
          }
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
