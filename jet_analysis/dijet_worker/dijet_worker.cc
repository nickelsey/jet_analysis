// dijet_worker.cc

#include "jet_analysis/dijet_worker/dijet_worker.hh"

#include "fastjet/tools/Subtractor.hh"

#include "TMath.h"

// default constructor
DijetWorker::DijetWorker() : DijetMatrix() { }

// single entry construction
DijetWorker::DijetWorker(fastjet::JetAlgorithm jet_alg_in,
                         double lead_pt_in,
                         double lead_R_in,
                         double sub_pt_in,
                         double sub_R_in,
                         double const_lead_pt_init_in,
                         double const_lead_pt_match_in,
                         double const_sub_pt_init_in,
                         double const_sub_pt_match_in,
                         double eta_in) :
      DijetMatrix(jet_alg_in, lead_pt_in, lead_R_in, sub_pt_in, sub_R_in,
                  const_lead_pt_init_in, const_lead_pt_match_in,
                  const_sub_pt_init_in, const_sub_pt_match_in, eta_in) { };

// set construction
DijetWorker::DijetWorker(std::set<fastjet::JetAlgorithm> jet_alg_in,
                         std::set<double> lead_pt_in,
                         std::set<double> lead_R_in,
                         std::set<double> sub_pt_in,
                         std::set<double> sub_R_in,
                         std::set<double> const_lead_pt_init_in,
                         std::set<double> const_lead_pt_match_in,
                         std::set<double> const_sub_pt_init_in,
                         std::set<double> const_sub_pt_match_in,
                         std::set<double> eta_in) :
      DijetMatrix(jet_alg_in, lead_pt_in, lead_R_in, sub_pt_in, sub_R_in,
                  const_lead_pt_init_in, const_lead_pt_match_in,
                  const_sub_pt_init_in, const_sub_pt_match_in, eta_in) { };

DijetWorker::DijetWorker(const DijetMatrix& rhs) :
      DijetMatrix(rhs) { }

std::unordered_map<std::string, ClusterOutput>& DijetWorker::Run(const std::vector<fastjet::PseudoJet>& input) {
  cluster_seq_lead_hard.clear();
  cluster_seq_sub_hard.clear();
  cluster_seq_lead_match.clear();
  cluster_seq_sub_match.clear();
  cluster_result.clear();
  
  // check to make sure the DijetMatrix has been initialized;
  // if not, do so
  if (Size() == 0)
    Initialize();
  std::cout << "NEW EVENT: " << std::endl;
  // loop over all dijet definitions and cluster;
  // there is only one sanity check for optimization currently -
  // if the hard or matching leading/subleading have the
  // same JetDefinition, it will not recluster for both leading &
  // subleading (possible 2x speedup). Further optimizations would be
  // possible comparing JetDefs across DijetDefinitions
  for (auto dijet_def : dijet_defs) {
    
    // get the key & leading/subleading definitions
    std::string key = dijet_def.first;
    std::shared_ptr<MatchDef> lead = dijet_def.second->lead;
    std::shared_ptr<MatchDef> sub = dijet_def.second->sub;
    std::cout << "KEY: " << key << std::endl;
    // make sure the dijet definition is valid
    if (!lead->IsValid() || !sub->IsValid())
      continue;
    
    ClusterOutput dijet_container(dijet_def.second);
    
    // do initial hard clustering for leading & subleading jets
    std::shared_ptr<fastjet::ClusterSequenceArea> cl_hard_lead = nullptr;
    std::shared_ptr<fastjet::ClusterSequenceArea> cl_hard_sub = nullptr;
    if (EquivalentClusterInput(lead->InitialJetDef(), sub->InitialJetDef())) {
      std::cout <<"constituent selector HARD: " << lead->InitialJetDef().ConstituentSelector().description() << std::endl;
      std::cout <<"jet def HARD: " << lead->InitialJetDef().description() << std::endl;
      std::cout <<"area def HARD: " << lead->InitialJetDef().AreaDefinition().description() << std::endl;
      cl_hard_lead = std::make_shared<fastjet::ClusterSequenceArea>(lead->InitialJetDef().ConstituentSelector()(input),
                                                                    lead->InitialJetDef(),
                                                                    lead->InitialJetDef().AreaDefinition());
      cl_hard_sub = cl_hard_lead;
    }
    else {
      cl_hard_lead = std::make_shared<fastjet::ClusterSequenceArea>(lead->InitialJetDef().ConstituentSelector()(input),
                                                                    lead->InitialJetDef(),
                                                                    lead->InitialJetDef().AreaDefinition());
      cl_hard_sub = std::make_shared<fastjet::ClusterSequenceArea>(sub->InitialJetDef().ConstituentSelector()(input),
                                                                   sub->InitialJetDef(),
                                                                   sub->InitialJetDef().AreaDefinition());
    }
    
    // insert into the dictionary
    cluster_seq_lead_hard.insert({key, cl_hard_lead});
    cluster_seq_sub_hard.insert({key, cl_hard_sub});
    
    // get output jets, and make sure neither are zero length
    std::vector<fastjet::PseudoJet> lead_hard_jets = fastjet::sorted_by_pt(lead->InitialJetDef().JetSelector()(cl_hard_lead->inclusive_jets()));
    std::vector<fastjet::PseudoJet> sublead_hard_jets = fastjet::sorted_by_pt(sub->InitialJetDef().JetSelector()(cl_hard_sub->inclusive_jets()));
    std::cout << "leading initial jet selector: " << lead->InitialJetDef().JetSelector().description() << std::endl;
    std::cout << "subleading initial jet selector: " << sub->InitialJetDef().JetSelector().description() << std::endl;
    std::cout << "leading jet count: " << lead_hard_jets.size() << std::endl;
    std::cout << "subleading jet count: " << sublead_hard_jets.size() << std::endl;
    
    if (lead_hard_jets.size() == 0)
      continue;
    
    // we have found a leading jet
    dijet_container.found_lead = true;
    
    // now, using the highest pT jet as the leading jet, try to find a
    // subleading jet that satisfies the dPhi cut
    fastjet::PseudoJet leading_hard_jet = lead_hard_jets[0];
    fastjet::Selector recoil_selector = !fastjet::SelectorRectangle(2.1, TMath::Pi() - dijet_def.second->dPhi);
    recoil_selector.set_reference(leading_hard_jet);
    std::vector<fastjet::PseudoJet> dphi_selected_recoil = fastjet::sorted_by_pt(recoil_selector(sublead_hard_jets));
    
    if (dphi_selected_recoil.size() == 0) {
      cluster_result.insert({key, dijet_container});
      continue;
    }
    
    // we have found a subleading jet that satisfies dPhi cuts
    dijet_container.found_sublead = true;
    
    fastjet::PseudoJet subleading_hard_jet = fastjet::sorted_by_pt(dphi_selected_recoil)[0];
    
    // now we will estimate the background for book-keeping
    // Energy density estimate from median ( pt_i / area_i )
    std::shared_ptr<fastjet::JetMedianBackgroundEstimator> lead_hard_bkg_est = nullptr;
    std::shared_ptr<fastjet::JetMedianBackgroundEstimator> sublead_hard_bkg_est = nullptr;
    if (EquivalentBkgEstimationInput(lead->InitialJetDef(), sub->InitialJetDef())) {
      lead_hard_bkg_est = std::make_shared<fastjet::JetMedianBackgroundEstimator>(lead->InitialJetDef().BackgroundSelector(),
                                                                                  lead->InitialJetDef().BackgroundJetDef(),
                                                                                  lead->InitialJetDef().BackgroundAreaDef());
      sublead_hard_bkg_est = lead_hard_bkg_est;
      lead_hard_bkg_est->set_particles(lead->InitialJetDef().ConstituentSelector()(input));
    }
    else {
      lead_hard_bkg_est = std::make_shared<fastjet::JetMedianBackgroundEstimator>(lead->InitialJetDef().BackgroundSelector(),
                                                                                  lead->InitialJetDef().BackgroundJetDef(),
                                                                                  lead->InitialJetDef().BackgroundAreaDef());
      sublead_hard_bkg_est = std::make_shared<fastjet::JetMedianBackgroundEstimator>(sub->InitialJetDef().BackgroundSelector(),
                                                                                     sub->InitialJetDef().BackgroundJetDef(),
                                                                                     sub->InitialJetDef().BackgroundAreaDef());
      lead_hard_bkg_est->set_particles(lead->InitialJetDef().ConstituentSelector()(input));
      sublead_hard_bkg_est->set_particles(sub->InitialJetDef().ConstituentSelector()(input));
    }
    
    // write the hard jets to the output container
    dijet_container.lead_hard = leading_hard_jet;
    dijet_container.lead_hard_rho = lead_hard_bkg_est->rho();
    dijet_container.lead_hard_sigma = lead_hard_bkg_est->sigma();
    dijet_container.sublead_hard = subleading_hard_jet;
    dijet_container.sublead_hard_rho = sublead_hard_bkg_est->rho();
    dijet_container.sublead_hard_sigma = sublead_hard_bkg_est->sigma();
    
    // now we will perform the matching jet finding
    // for both leading & subleading jets
    std::shared_ptr<fastjet::ClusterSequenceArea> cl_match_lead = nullptr;
    std::shared_ptr<fastjet::ClusterSequenceArea> cl_match_sub = nullptr;
    
    if (EquivalentClusterInput(lead->MatchedJetDef(), sub->MatchedJetDef())) {
      std::cout <<"constituent selector MATCH: " << lead->MatchedJetDef().ConstituentSelector().description() << std::endl;
      std::cout <<"jet def MATCH: " << lead->MatchedJetDef().description() << std::endl;
      std::cout <<"area def MATCH: " << lead->MatchedJetDef().AreaDefinition().description() << std::endl;
      cl_match_lead = std::make_shared<fastjet::ClusterSequenceArea>(lead->MatchedJetDef().ConstituentSelector()(input),
                                                                     lead->MatchedJetDef(),
                                                                     lead->MatchedJetDef().AreaDefinition());
      cl_match_sub = cl_match_lead;
    }
    else {
      cl_match_lead = std::make_shared<fastjet::ClusterSequenceArea>(lead->MatchedJetDef().ConstituentSelector()(input),
                                                                     lead->MatchedJetDef(),
                                                                     lead->MatchedJetDef().AreaDefinition());
      cl_match_sub = std::make_shared<fastjet::ClusterSequenceArea>(sub->MatchedJetDef().ConstituentSelector()(input),
                                                                    sub->MatchedJetDef(),
                                                                    sub->MatchedJetDef().AreaDefinition());
    }
    
    // insert into the dictionary
    cluster_seq_lead_match.insert({key, cl_match_lead});
    cluster_seq_sub_match.insert({key, cl_match_sub});
    
    // get the resulting jets and make sure neither are zero length
    std::vector<fastjet::PseudoJet> lead_match_jets = fastjet::sorted_by_pt(lead->MatchedJetDef().JetSelector()(cl_match_lead->inclusive_jets()));
    std::vector<fastjet::PseudoJet> sublead_match_jets = fastjet::sorted_by_pt(sub->MatchedJetDef().JetSelector()(cl_match_sub->inclusive_jets()));
    std::cout << "leading matched jet selector: " << lead->MatchedJetDef().JetSelector().description() << std::endl;
    std::cout << "subleading matched jet selector: " << sub->MatchedJetDef().JetSelector().description() << std::endl;
    std::cout << "leading matched jet count: " << lead_match_jets.size() << std::endl;
    std::cout << "subleading matched jet count: " << sublead_match_jets.size() << std::endl;
    
    if (lead_match_jets.size() == 0 || sublead_match_jets.size() == 0) {
      cluster_result.insert({key, dijet_container});
      continue;
    }
    
    // now do background subtraction
    std::shared_ptr<fastjet::JetMedianBackgroundEstimator> lead_match_bkg_est = nullptr;
    std::shared_ptr<fastjet::JetMedianBackgroundEstimator> sublead_match_bkg_est = nullptr;
    std::vector<fastjet::PseudoJet> lead_subtracted_jets;
    std::vector<fastjet::PseudoJet> sublead_subtracted_jets;
    
    if (EquivalentBkgEstimationInput(lead->MatchedJetDef(), sub->MatchedJetDef())) {
      lead_match_bkg_est = std::make_shared<fastjet::JetMedianBackgroundEstimator>(lead->MatchedJetDef().BackgroundSelector(),
                                                                                   lead->MatchedJetDef().BackgroundJetDef(),
                                                                                   lead->MatchedJetDef().BackgroundAreaDef());
      sublead_match_bkg_est = lead_match_bkg_est;
      lead_match_bkg_est->set_particles(lead->MatchedJetDef().ConstituentSelector()(input));
      fastjet::Subtractor bkgdSubtractor(&(*lead_match_bkg_est));
      lead_subtracted_jets = fastjet::sorted_by_pt(bkgdSubtractor(cl_match_lead->inclusive_jets()));
      sublead_subtracted_jets = fastjet::sorted_by_pt(bkgdSubtractor(cl_match_sub->inclusive_jets()));
    }
    else {
      lead_match_bkg_est = std::make_shared<fastjet::JetMedianBackgroundEstimator>(lead->MatchedJetDef().BackgroundSelector(),
                                                                                  lead->MatchedJetDef().BackgroundJetDef(),
                                                                                  lead->MatchedJetDef().BackgroundAreaDef());
      sublead_match_bkg_est = std::make_shared<fastjet::JetMedianBackgroundEstimator>(sub->MatchedJetDef().BackgroundSelector(),
                                                                                     sub->MatchedJetDef().BackgroundJetDef(),
                                                                                     sub->MatchedJetDef().BackgroundAreaDef());
      lead_match_bkg_est->set_particles(lead->MatchedJetDef().ConstituentSelector()(input));
      sublead_match_bkg_est->set_particles(sub->MatchedJetDef().ConstituentSelector()(input));
      fastjet::Subtractor lead_subtractor(&(*lead_match_bkg_est));
      fastjet::Subtractor sublead_subtractor(&(*sublead_match_bkg_est));
      lead_subtracted_jets = fastjet::sorted_by_pt(lead_subtractor(cl_match_lead->inclusive_jets()));
      sublead_subtracted_jets = fastjet::sorted_by_pt(sublead_subtractor(cl_match_sub->inclusive_jets()));
    }
  
    // now match to the hard jets
    fastjet::Selector lead_circle_selector = fastjet::SelectorCircle(lead->dR());
    lead_circle_selector.set_reference(leading_hard_jet);
    fastjet::Selector sublead_circle_selector = fastjet::SelectorCircle(sub->dR());
    sublead_circle_selector.set_reference(subleading_hard_jet);
    
    std::vector<fastjet::PseudoJet> matched_to_leading = fastjet::sorted_by_pt(lead_circle_selector(lead_subtracted_jets));
    std::vector<fastjet::PseudoJet> matched_to_subleading = fastjet::sorted_by_pt(sublead_circle_selector(sublead_subtracted_jets));
    
    
    // make sure there is a valid matched jet
    if (matched_to_leading.size() == 0 || matched_to_subleading.size() == 0) {
      cluster_result.insert({key, dijet_container});
      continue;
    }
    
    dijet_container.found_match = true;
    
    fastjet::PseudoJet leading_matched_jet = matched_to_leading[0];
    fastjet::PseudoJet subleading_matched_jet = matched_to_subleading[0];
    
    // put into output container
    dijet_container.lead_match = leading_matched_jet;
    dijet_container.lead_match_rho = lead_match_bkg_est->rho();
    dijet_container.lead_match_sigma = lead_match_bkg_est->sigma();
    
    dijet_container.sublead_match = subleading_matched_jet;
    dijet_container.sublead_match_rho = sublead_match_bkg_est->rho();
    dijet_container.sublead_match_sigma = sublead_match_bkg_est->sigma();
    
    // and put the container in the output
    cluster_result.insert({key, dijet_container});
    
  }
  
  return cluster_result;
}

bool DijetWorker::EquivalentAreaDefinition(const fastjet::AreaDefinition& a1, const fastjet::AreaDefinition& a2) {
  
  // compare area definitions
  if (a1.area_type() != a2.area_type() ||
      a1.ghost_spec().ghost_maxrap() != a2.ghost_spec().ghost_maxrap() ||
      a1.ghost_spec().repeat() != a2.ghost_spec().repeat() ||
      a1.ghost_spec().ghost_area() != a2.ghost_spec().ghost_area() ||
      a1.ghost_spec().grid_scatter() != a2.ghost_spec().grid_scatter() ||
      a1.ghost_spec().pt_scatter() != a2.ghost_spec().pt_scatter() ||
      a1.ghost_spec().mean_ghost_pt() != a2.ghost_spec().mean_ghost_pt())
    return false;
  
  return true;
}

bool DijetWorker::EquivalentClusterInput(const JetDef& c1, const JetDef& c2) {
  // compare the jet definition
  if (c1.R() != c2.R() ||
      c1.jet_algorithm() != c2.jet_algorithm() ||
      c1.strategy() != c2.strategy() ||
      c1.recombination_scheme() != c2.recombination_scheme())
    return false;
  
  // compare AreaDefinitions
  if (!EquivalentAreaDefinition(c1.AreaDefinition(), c2.AreaDefinition()))
    return false;
  
  // compare constituent selector since I dont want to
  // compare all possible selectors, will only compare
  // description strings, so if they werent constructed
  // identically this will fail, even if they give
  // identical output
  if (c1.ConstituentSelector().description() != c2.ConstituentSelector().description())
    return false;
  
  return true;
}

bool DijetWorker::EquivalentBkgEstimationInput(const JetDef& c1, const JetDef& c2) {
  
  // first. compare constituent and bkg jet selector
  if (c1.ConstituentSelector().description() != c2.ConstituentSelector().description())
    return false;
  
  if (c1.BackgroundSelector().description() != c2.BackgroundSelector().description())
    return false;
  
  // compare background area definitions
  if (!EquivalentAreaDefinition(c1.BackgroundAreaDef(), c2.BackgroundAreaDef()))
    return false;
  
  // finally compare background jet definitions
  if (c1.BackgroundJetDef().R() != c2.BackgroundJetDef().R() ||
      c1.BackgroundJetDef().jet_algorithm() != c2.BackgroundJetDef().jet_algorithm() ||
      c1.BackgroundJetDef().recombination_scheme() != c2.BackgroundJetDef().recombination_scheme() ||
      c1.BackgroundJetDef().strategy() != c2.BackgroundJetDef().strategy())
    return false;
  
  return true;
}
