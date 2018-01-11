// differential_aj.cxx

#include <iostream>
#include <string>
#include <fstream>
#include <unordered_map>
#include <random>
#include <exception>

#include "jet_analysis/util/arg_helper.hh"
#include "jet_analysis/util/trigger_lookup.hh"
#include "jet_analysis/util/reader_util.hh"
#include "jet_analysis/util/vector_conversion.hh"
#include "jet_analysis/efficiency/run4_eff.hh"
#include "jet_analysis/dijet_worker/dijet_worker.hh"

#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TH2.h"
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
  string name      = "job"; /* output file name */
  string id        = "0";   /* job id */
  string input     = "";    /* root file/root file list*/
  string embed     = "";    /* root file/file list for embedding*/
  int    reuse     = 1;     /* number of events an HT event should be embedded in */
  string reader    = "";    /* settings file for primary reader */
  string readeremb = "";    /* settings for embedding reader */
  string out_dir   = "";    /* directory to save output in */
  string tow_list  = "";    /* list of hot towers to remove */
  string run_list  = "";    /* list of runs to remove */
  string triggers  = "";    /* triggers to consider (see trigger_lookup.hh) */
  string trigemb   = "";    /* triggers to consider for embedding (should be mb) */
  string const_eta = "";    /* constituent eta cut */
  string lc_pt     = "";    /* leading hard constituent pt cut */
  string sc_pt     = "";    /* subleading hard constituent pt cut */
  string lcm_pt    = "";    /* leading matched constituent pt cut */
  string scm_pt    = "";    /* subleading matched constituent pt cut */
  string lead_r    = "";    /* lead jet radii */
  string sub_r     = "";    /* sublead jet radii */
  string lj_pt     = "";    /* leading hard jet pt cut */
  string sj_pt     = "";    /* subleading hard jet pt cut */
  int apply_effic  =  0;    /* whether or not to do efficiency corrections */
};

int main(int argc, char* argv[]) {
  
  // parse command line options
  Options opts;
  for (int i = 1; i < argc; ++i) {
    if (ParseStrFlag(string(argv[i]), "--name", &opts.name) ||
        ParseStrFlag(string(argv[i]), "--id", &opts.id) ||
        ParseStrFlag(string(argv[i]), "--input", &opts.input) ||
        ParseStrFlag(string(argv[i]), "--embed", &opts.embed) ||
        ParseIntFlag(string(argv[i]), "--reuseTrigger", &opts.reuse) ||
        ParseStrFlag(string(argv[i]), "--readerSetting", &opts.reader) ||
        ParseStrFlag(string(argv[i]), "--embedReaderSetting", &opts.readeremb) ||
        ParseStrFlag(string(argv[i]), "--outDir", &opts.out_dir) ||
        ParseStrFlag(string(argv[i]), "--towList", &opts.tow_list) ||
        ParseStrFlag(string(argv[i]), "--runList", &opts.run_list) ||
        ParseStrFlag(string(argv[i]), "--triggers", &opts.triggers) ||
        ParseStrFlag(string(argv[i]), "--embedTriggers", &opts.trigemb) ||
        ParseIntFlag(string(argv[i]), "--efficiency", &opts.apply_effic) ||
        ParseStrFlag(string(argv[i]), "--constEta", &opts.const_eta) ||
        ParseStrFlag(string(argv[i]), "--leadConstPt", &opts.lc_pt) ||
        ParseStrFlag(string(argv[i]), "--subConstPt", &opts.sc_pt) ||
        ParseStrFlag(string(argv[i]), "--leadConstPtMatch", &opts.lcm_pt) ||
        ParseStrFlag(string(argv[i]), "--subConstPtMatch", &opts.scm_pt) ||
        ParseStrFlag(string(argv[i]), "--leadR", &opts.lead_r) ||
        ParseStrFlag(string(argv[i]), "--subR", &opts.sub_r) ||
        ParseStrFlag(string(argv[i]), "--leadJetPt", &opts.lj_pt) ||
        ParseStrFlag(string(argv[i]), "--subJetPt", &opts.sj_pt)) continue;
    std::cerr << "Unknown command line option: " << argv[i] << std::endl;
    return 1;
  }
  
  // sanity check - if no embedding is specified,
  // triggers can't be reused (just double counting at
  // that point). Also make sure its >= 1, so it doesn't
  // break...
  if (opts.embed.empty() || opts.reuse < 1)
    opts.reuse = 1;
  
  // check to make sure the input file paths are sane
  if (!boost::filesystem::exists(opts.input)) {
    std::cerr << "input root file does not exist: ";
    std::cerr << opts.input << std::endl;
    std::cerr << "exiting" << std::endl;
    return 1;
  }
  
  if (!boost::filesystem::exists(opts.embed) && !opts.embed.empty()) {
    std::cerr << "input embedding file specified but doesn't exist:";
    std::cerr << opts.embed << std::endl;
    std::cerr << "exiting" << std::endl;
    return 1;
  }
  
  // first, build our input chains -
  // the default chain must be initialized, but
  // no embedding is necessary
  TChain* chain = NewChainFromInput(opts.input);
  TChain* embed_chain = NewChainFromInput(opts.embed);
  
  // build output directory if it doesn't exist
  // using the boost::filesystem library
  if (opts.out_dir.empty())
    opts.out_dir = "tmp";
  boost::filesystem::path dir(opts.out_dir.c_str());
  boost::filesystem::create_directories(dir);
  
  // create output file from the given directory, name & id
  string outfile_name = opts.out_dir + "/" + opts.name + opts.id + ".root";
  TFile out(outfile_name.c_str(), "RECREATE");
  
  // initialize the reader(s)
  TStarJetPicoReader* reader = new TStarJetPicoReader();
  if (!opts.reader.empty()) {
    InitReader(reader, chain, opts.reader, opts.tow_list, opts.run_list);
  }
  else {
    InitReaderWithDefaults(reader, chain, opts.tow_list, opts.run_list);
  }
  
  // initialize the embedding reader if an embedding
  // chain has been created
  TStarJetPicoReader* reader_embed = nullptr;
  if (embed_chain != nullptr) {
    reader_embed = new TStarJetPicoReader();
    if (!opts.readeremb.empty()) {
      InitReader(reader_embed, embed_chain, opts.readeremb,
                 opts.tow_list, opts.run_list);
    }
    else {
      InitReaderWithDefaults(reader_embed, embed_chain, opts.tow_list,
                             opts.run_list);
    }
    
    // I currently use all MB available data for embedding, so we want to start
    // each job with a different initial embedding event
    double prob_min = 0.0;
    double prob_max = reader_embed->GetNOfEvents() - 1;
    std::random_device generator;
    std::uniform_int_distribution<> flat_probability(prob_min, prob_max);
    reader_embed->ReadEvent(flat_probability(generator));
  }
  
  // get the triggers IDs that will be used
  std::set<unsigned> triggers = GetTriggerIDs(opts.triggers);
  std::set<unsigned> embed_triggers;
  if (!opts.trigemb.empty())
    embed_triggers = GetTriggerIDs(opts.trigemb);
  
  std::cout << "taking triggers: " << opts.triggers << " for primary" << std::endl;
  std::cout << "trigger ids: ";
  for (auto i : triggers)
    std::cout << i << " ";
  std::cout << std::endl;
  
  if (reader_embed != nullptr) {
    std::cout << "taking triggers: " << opts.trigemb << " for embedding" << std::endl;
    std::cout << "trigger ids: ";
    for (auto i : embed_triggers)
      std::cout << i << " ";
    std::cout << std::endl;
  }
  
  // now parse analysis variables from command line arguments
  
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
  std::cout << "worker initialized - number of dijet definitions: " << worker.Size() <<std::endl;
  std::set<std::string> keys = worker.Keys();
  for (auto key : keys)
    std::cout << key << std::endl;
  
  // create an output tree for each definition
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
    
    if (reader_embed != nullptr) {
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
    
    if (reader_embed != nullptr) {
      tmp->Branch("embed_eventid", &embed_eventid_dict[key]);
      tmp->Branch("embed_runid", &embed_runid_dict[key]);
      tmp->Branch("embed_refmult", &embed_refmult_dict[key]);
      tmp->Branch("embed_grefmult", &embed_grefmult_dict[key]);
      tmp->Branch("embed_refmultcorr", &embed_refmultcorr_dict[key]);
      tmp->Branch("embed_grefmultcorr", &embed_grefmultcorr_dict[key]);
      tmp->Branch("embed_cent", &embed_cent_dict[key]);
      tmp->Branch("embed_rp", &embed_rp_dict[key]);
      tmp->Branch("embed_zdcrate", &embed_zdc_dict[key]);
      tmp->Branch("embed_vz", &embed_vz_dict[key]);
    }
    
    trees.insert({key, tmp});
  }
  
  // define histograms for coincidence measurements
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
  
  // initialize the efficiency class
  Run4Eff efficiency;
  
  // for now, use hard coded centrality definition for run 14
  std::vector<int> refcent_def{420, 364, 276, 212, 156, 108, 68, 44, 28, 12, 0};
  
  // when getting ratios of efficiencies, we need to specify the comparison bin
  // not offered as a command line option, because I don't see it changing from
  // the most central bin anytime soon
  const int refcent_reference = 0;
  
  // random number distribution for a flat probability
  // used when discarding tracks based on relative efficiency
  double prob_min = 0.0;
  double prob_max = 1.0;
  std::uniform_real_distribution<> flat_probability(prob_min, prob_max);
  // create the generator
  // standard mersenne_twister_engine seeded with a constant,
  // so that it is reproducible
  std::mt19937 gen(14342);
  
  // define a selector to reject low momentum tracks
  fastjet::Selector track_pt_min_selector = fastjet::SelectorPtMin(0.2) && fastjet::SelectorAbsRapMax(1.0);
  
  // start the analysis loop
  try {
    while (reader->NextEvent()) {
      // Print out reader status every 10 seconds
      reader->PrintStatus(10);
      
      // headers for convenience
      TStarJetPicoEventHeader* header = reader->GetEvent()->GetHeader();
      TStarJetPicoEventHeader* header_embed = nullptr;
      if (reader_embed != nullptr) {
        header_embed = reader_embed->GetEvent()->GetHeader();
      }
      
      // check if event fired a trigger we will use
      if (triggers.size() != 0) {
        bool use_event = false;
        for (auto trigger : triggers)
          if (header->HasTriggerId(trigger))
            use_event = true;
        if (!use_event) continue;
      }
      
      // get event reference centrality
      int centrality = 0;
      int eff_corr_cent = 0; // used for efficiency corrections
      for (unsigned i = 0; i < refcent_def.size(); ++i)
        if (header->GetReferenceMultiplicity() > refcent_def[i]) {
          centrality = i;
          break;
        }
      if (centrality > 8)
        eff_corr_cent = 8;
      else
        eff_corr_cent = centrality;
      
      // get the vector container
      TStarJetVectorContainer<TStarJetVector>* container = reader->GetOutputContainer();
      std::vector<fastjet::PseudoJet> primary_particles;
      ConvertTStarJetVector(container, primary_particles);
      
      // select tracks above the minimum pt threshold
      primary_particles = track_pt_min_selector(primary_particles);
      
      // now, we will loop over the same data opts.reuse times, using the same
      // triggered event. This will give us the ability to generate multiple
      // embedding events per triggered event, giving higher statistics for our
      // pp reference
      
      for (int i = 0; i < opts.reuse; ++i) {
        
        // container that will be passed to the DijetWorker
        std::vector<fastjet::PseudoJet> input;
        
        // if relative tracking efficiency smearing is requested,
        // discard tracks randomly by the ratio of efficiencies
        // either between AuAu of different centralities, or between
        // AuAu of a specific centrality and PP
        if (opts.apply_effic) {
          if (reader_embed == nullptr) {
            for (auto vec : primary_particles) {
              if (efficiency.CentRatio(vec.pt(), vec.eta(), eff_corr_cent, refcent_reference) > flat_probability(gen))
                  input.push_back(vec);
            }
          }
          else {
            for (auto vec : primary_particles) {
              if (efficiency.AuAuPPRatio(vec.pt(), vec.eta(), refcent_reference) > flat_probability(gen))
                input.push_back(vec);
            }
          }
        }
        // if no efficiency correction is being applied, copy the full set of particles
        else {
          input = std::vector<fastjet::PseudoJet>(primary_particles.begin(),
                                                  primary_particles.end());
        }
        
        int embed_centrality = 0;
        int eff_corr_embed_cent = 0;
        // if embedding is requested, add to the input
        if (reader_embed != nullptr) {
          // read next event, loop to event 0 if needed
          if (!GetNextValidEvent(reader_embed, embed_triggers)) {
            std::cerr << "no events found for embedding, given trigger requirements: exiting" << std::endl;
            return 1;
          }
        
        // get the centrality definition of the embedding data
        for (unsigned i = 0; i < refcent_def.size(); ++i)
          if (header_embed->GetReferenceMultiplicity() > refcent_def[i]) {
            embed_centrality = i;
            break;
          }
          if (embed_centrality > 8)
            eff_corr_embed_cent = 8;
          else
            eff_corr_embed_cent = embed_centrality;
        
          // get the vector container for embedding
          TStarJetVectorContainer<TStarJetVector>* container_embed = reader_embed->GetOutputContainer();
          std::vector<fastjet::PseudoJet> embed_particles;
          ConvertTStarJetVector(container_embed, embed_particles);
        
          // again, apply minimal kinematic cuts
          embed_particles = track_pt_min_selector(embed_particles);
        
          // if relative tracking efficiency smearing is requested,
          // discard tracks randomly by the ratio of efficiencies
          // it is assumed that the embedding event is AuAu
          if (opts.apply_effic) {
            for (auto vec : embed_particles)
              if (efficiency.CentRatio(vec.pt(), vec.eta(), eff_corr_embed_cent, refcent_reference) > flat_probability(gen))
                input.push_back(vec);
          }
          else {
            input.insert(input.end(), embed_particles.begin(), embed_particles.end());
          }
        }
        
        // run the worker
        auto worker_out = worker.Run(input);
        
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
            refmultcorr_dict[key] = header->GetCorrectedReferenceMultiplicity();
            grefmultcorr_dict[key] = header->GetCorrectedGReferenceMultiplicity();
            cent_dict[key] = centrality;
            zdcrate_dict[key] = header->GetZdcCoincidenceRate();
            reactionplane_dict[key] = header->GetReactionPlaneAngle();
            nglobal_dict[key] = header->GetNGlobalTracks();

            // if embedding is being done
            if (reader_embed != nullptr) {
              embed_runid_dict[key] = header_embed->GetRunId();
              embed_eventid_dict[key] = header_embed->GetEventId();
              embed_refmult_dict[key] = header_embed->GetReferenceMultiplicity();
              embed_grefmult_dict[key] = header_embed->GetGReferenceMultiplicity();
              embed_refmultcorr_dict[key] = header_embed->GetCorrectedReferenceMultiplicity();
              embed_grefmultcorr_dict[key] = header_embed->GetCorrectedGReferenceMultiplicity();
              embed_cent_dict[key] = embed_centrality;
              embed_vz_dict[key] = header_embed->GetPrimaryVertexZ();
              embed_zdc_dict[key] = header_embed->GetZdcCoincidenceRate();
              embed_rp_dict[key] = header_embed->GetReactionPlaneAngle();
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
