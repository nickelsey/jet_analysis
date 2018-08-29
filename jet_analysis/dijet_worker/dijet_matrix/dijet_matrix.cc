// dijet_matrix.cc

#include "jet_analysis/dijet_worker/dijet_matrix/dijet_matrix.hh"

#include "jet_analysis/util/selector_compare.hh"
#include "jet_analysis/util/dijet_key.hh"

// default constructor
DijetMatrix::DijetMatrix() :
      force_constituent_pt_equality(true),
      force_constituent_eta_equality(true),
      strategy(fastjet::Best),
      scheme(fastjet::E_scheme),
      area_type(fastjet::active_area_explicit_ghosts),
      ghost_repeat(fastjet::gas::def_repeat),
      ghost_area(fastjet::gas::def_ghost_area),
      grid_scatter(fastjet::gas::def_grid_scatter),
      pt_scatter(fastjet::gas::def_pt_scatter),
      mean_ghost_pt(fastjet::gas::def_mean_ghost_pt),
      bkg_definition(fastjet::JetDefinition(fastjet::kt_algorithm, 0.4)) { }

// single entry construction
DijetMatrix::DijetMatrix(fastjet::JetAlgorithm jet_alg_in,
                         double lead_pt_in,
                         double lead_R_in,
                         double sub_pt_in,
                         double sub_R_in,
                         double const_lead_pt_init_in,
                         double const_lead_pt_match_in,
                         double const_sub_pt_init_in,
                         double const_sub_pt_match_in,
                         double eta_in) :
      const_eta(std::set<double>{eta_in}),
      const_lead_pt_init(std::set<double>{const_lead_pt_init_in}),
      const_lead_pt_match(std::set<double>{const_lead_pt_match_in}),
      const_sub_pt_init(std::set<double>{const_sub_pt_init_in}),
      const_sub_pt_match(std::set<double>{const_sub_pt_match_in}),
      jet_algorithm(std::set<fastjet::JetAlgorithm>{jet_alg_in}),
      lead_pt(std::set<double>{lead_pt_in}),
      lead_R(std::set<double>{lead_R_in}),
      sub_pt(std::set<double>{sub_pt_in}),
      sub_R(std::set<double>{sub_R_in}),
      force_constituent_pt_equality(true),
      force_constituent_eta_equality(true),
      strategy(fastjet::Best),
      scheme(fastjet::E_scheme),
      area_type(fastjet::active_area_explicit_ghosts),
      ghost_repeat(fastjet::gas::def_repeat),
      ghost_area(fastjet::gas::def_ghost_area),
      grid_scatter(fastjet::gas::def_grid_scatter),
      pt_scatter(fastjet::gas::def_pt_scatter),
      mean_ghost_pt(fastjet::gas::def_mean_ghost_pt),
      bkg_definition(fastjet::JetDefinition(fastjet::kt_algorithm, 0.4))
      { }


// set construction
DijetMatrix::DijetMatrix(std::set<fastjet::JetAlgorithm> jet_alg_in,
                         std::set<double> lead_pt_in,
                         std::set<double> lead_R_in,
                         std::set<double> sub_pt_in,
                         std::set<double> sub_R_in,
                         std::set<double> const_lead_pt_init_in,
                         std::set<double> const_lead_pt_match_in,
                         std::set<double> const_sub_pt_init_in,
                         std::set<double> const_sub_pt_match_in,
                         std::set<double> eta_in) :
      const_eta(eta_in),
      const_lead_pt_init(const_lead_pt_init_in),
      const_lead_pt_match(const_lead_pt_match_in),
      const_sub_pt_init(const_sub_pt_init_in),
      const_sub_pt_match(const_sub_pt_match_in),
      jet_algorithm(jet_alg_in),
      lead_pt(lead_pt_in),
      lead_R(lead_R_in),
      sub_pt(sub_pt_in),
      sub_R(sub_R_in),
      force_constituent_pt_equality(true),
      force_constituent_eta_equality(true),
      strategy(fastjet::Best),
      scheme(fastjet::E_scheme),
      area_type(fastjet::active_area_explicit_ghosts),
      ghost_repeat(fastjet::gas::def_repeat),
      ghost_area(fastjet::gas::def_ghost_area),
      grid_scatter(fastjet::gas::def_grid_scatter),
      pt_scatter(fastjet::gas::def_pt_scatter),
      mean_ghost_pt(fastjet::gas::def_mean_ghost_pt),
      bkg_definition(fastjet::JetDefinition(fastjet::kt_algorithm, 0.4))
      { }

DijetMatrix::DijetMatrix(const DijetMatrix& rhs) :
      const_eta(rhs.const_eta),
      const_lead_pt_init(rhs.const_lead_pt_init),
      const_lead_pt_match(rhs.const_lead_pt_match),
      const_sub_pt_init(rhs.const_sub_pt_init),
      const_sub_pt_match(rhs.const_sub_pt_init),
      jet_algorithm(rhs.jet_algorithm),
      lead_pt(rhs.lead_pt),
      lead_R(rhs.lead_R),
      sub_pt(rhs.sub_pt),
      sub_R(rhs.sub_R),
      force_constituent_pt_equality(rhs.force_constituent_pt_equality),
      force_constituent_eta_equality(rhs.force_constituent_eta_equality),
      strategy(rhs.strategy),
      scheme(rhs.scheme),
      area_type(rhs.area_type),
      ghost_repeat(rhs.ghost_repeat),
      ghost_area(rhs.ghost_area),
      grid_scatter(rhs.grid_scatter),
      pt_scatter(rhs.pt_scatter),
      mean_ghost_pt(rhs.mean_ghost_pt),
      bkg_definition(rhs.bkg_definition)
      { }

void DijetMatrix::Clear() {
  keys.clear();
  const_eta.clear();
  const_lead_pt_init.clear();
  const_sub_pt_init.clear();
  const_lead_pt_match.clear();
  const_sub_pt_match.clear();
  dijet_defs.clear();
  lead_matchdefs.clear();
  sub_matchdefs.clear();
  jet_algorithm.clear();
  lead_pt.clear();
  lead_R.clear();
  sub_pt.clear();
  sub_R.clear();
}

void DijetMatrix::ClearDijetDefs() {
  dijet_defs.clear();
  lead_matchdefs.clear();
  sub_matchdefs.clear();
  keys.clear();
}

void DijetMatrix::ForceConstituentPtEquality(bool flag) {
  force_constituent_pt_equality = flag;
  CheckToUpdate();
}

void DijetMatrix::ForceConstituentEtaEquality(bool flag) {
  force_constituent_eta_equality = flag;
  CheckToUpdate();
}

void DijetMatrix::AddConstituentEta(double eta) {
  const_eta.insert(eta);
  CheckToUpdate();
}
void DijetMatrix::AddConstituentEta(std::set<double> eta) {
  for (auto& val : eta) {
    const_eta.insert(val);
  }
  CheckToUpdate();
}

void DijetMatrix::AddConstituentLeadInitialPt(double pt) {
  const_lead_pt_init.insert(pt);
  CheckToUpdate();
}

void DijetMatrix::AddConstituentLeadInitialPt(std::set<double> pt) {
  for (auto& val : pt) {
    const_lead_pt_init.insert(val);
  }
  CheckToUpdate();
}

void DijetMatrix::AddConstituentLeadMatchPt(double pt) {
  const_lead_pt_match.insert(pt);
  CheckToUpdate();
}

void DijetMatrix::AddConstituentLeadMatchPt(std::set<double> pt) {
  for (auto& val : pt) {
    const_lead_pt_match.insert(val);
  }
  CheckToUpdate();
}

void DijetMatrix::AddConstituentSubInitialPt(double pt) {
  const_sub_pt_init.insert(pt);
  CheckToUpdate();
}

void DijetMatrix::AddConstituentSubInitialPt(std::set<double> pt) {
  for (auto& val : pt) {
    const_sub_pt_init.insert(val);
  }
  CheckToUpdate();
}

void DijetMatrix::AddConstituentSubMatchPt(double pt) {
  const_sub_pt_match.insert(pt);
  CheckToUpdate();
}

void DijetMatrix::AddConstituentSubMatchPt(std::set<double> pt) {
  for (auto& val : pt) {
    const_sub_pt_match.insert(val);
  }
  CheckToUpdate();
}

void DijetMatrix::AddJetAlgorithm(fastjet::JetAlgorithm alg) {
  jet_algorithm.insert(alg);
  CheckToUpdate();
}

void DijetMatrix::AddJetAlgorithm(std::set<fastjet::JetAlgorithm> alg) {
  for (auto& val : alg) {
    jet_algorithm.insert(val);
  }
  CheckToUpdate();
}

void DijetMatrix::AddLeadJetPt(double pt) {
  lead_pt.insert(pt);
  CheckToUpdate();
}

void DijetMatrix::AddLeadJetPt(std::set<double> pt) {
  for (auto& val : pt) {
    lead_pt.insert(val);
  }
  CheckToUpdate();
}

void DijetMatrix::AddLeadJetR(double R) {
  lead_R.insert(R);
  CheckToUpdate();
}

void DijetMatrix::AddLeadJetR(std::set<double> R) {
  for (auto& val : R) {
    lead_R.insert(val);
  }
  CheckToUpdate();
}

void DijetMatrix::AddSubJetPt(double pt) {
  sub_pt.insert(pt);
  CheckToUpdate();
}

void DijetMatrix::AddSubJetPt(std::set<double> pt) {
  for (auto& val : pt) {
    sub_pt.insert(val);
  }
  CheckToUpdate();
}

void DijetMatrix::AddSubJetR(double R) {
  sub_R.insert(R);
  CheckToUpdate();
}

void DijetMatrix::AddSubJetR(std::set<double> R) {
  for (auto& val : R) {
    sub_R.insert(val);
  }
  CheckToUpdate();
}

void DijetMatrix::Initialize() {
  // if already initialzed, remove old definitions
  if (dijet_defs.size() != 0) {
    ClearDijetDefs();
  }
  
  // check to make sure there is at least one valid dijet
  // definition, either through custom definitions,
  // or parameters. Create default parameters otherwise
  InitializeEmptyFields();
  
  // first step - get our leading & subleading jet definitions
  std::vector<fastjet::JetDefinition> lead_jet_defs = FillLeadJetDefinitions();
  std::vector<fastjet::JetDefinition> sub_jet_defs = FillSubJetDefinitions();
  
  // start creating the leading jet MatchDefs
  for (auto& lead_fj_def : lead_jet_defs) {
    for (auto& lead_jet_pt : lead_pt) {
      for (auto& lead_const_pt_init : const_lead_pt_init) {
        for (auto& lead_const_pt_match : const_lead_pt_match) {
          for (auto& lead_const_eta : const_eta) {
  
            // get a few needed parameters
            double R = lead_fj_def.R();
            double bkg_R = bkg_definition.R();
            double jet_eta_max = lead_const_eta - R;
            double bkg_jet_eta_max = lead_const_eta - bkg_R;
            
            // build constituent & jet selectors
            fastjet::Selector init_const_selector = fastjet::SelectorAbsRapMax(lead_const_eta)
                                                  && fastjet::SelectorPtMin(lead_const_pt_init);
            fastjet::Selector match_const_selector = fastjet::SelectorAbsRapMax(lead_const_eta)
                                                   && fastjet::SelectorPtMin(lead_const_pt_match);
            fastjet::Selector init_jet_selector = fastjet::SelectorAbsRapMax(jet_eta_max)
                                                && fastjet::SelectorPtMin(lead_jet_pt);
            fastjet::Selector match_jet_selector = fastjet::SelectorIdentity();
              
            fastjet::GhostedAreaSpec ghost_def(jet_eta_max + 2 * R, ghost_repeat, ghost_area,
                                               grid_scatter, pt_scatter, mean_ghost_pt);
            fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts,
                                              ghost_def);
            
            fastjet::Selector bkg_selector = fastjet::SelectorAbsRapMax(lead_const_eta - bkg_R)
                                           && (!fastjet::SelectorNHardest(2));
            fastjet::GhostedAreaSpec bkg_ghost_def(bkg_jet_eta_max + 2 * bkg_R, ghost_repeat, ghost_area,
                                                   grid_scatter, pt_scatter, mean_ghost_pt);
            fastjet::AreaDefinition bkg_area_def(fastjet::active_area_explicit_ghosts,
                                                 bkg_ghost_def);
            // create the initial JetDef & the
            // matched JetDef def
            JetDef init_def(lead_fj_def, area_def, bkg_definition, bkg_area_def, bkg_selector);
            init_def.SetConstituentSelector(init_const_selector);
            init_def.SetJetSelector(init_jet_selector);
            JetDef match_def(lead_fj_def, area_def, bkg_definition, bkg_area_def, bkg_selector);
            match_def.SetConstituentSelector(match_const_selector);
            match_def.SetJetSelector(match_jet_selector);
              
            // create the MatchDef
            lead_matchdefs.insert(std::make_unique<MatchDef>(init_def, match_def));
          }
        }
      }
    }
  }
  
  // do the same for subleading
  for (auto& sub_fj_def : sub_jet_defs) {
    for (auto& sub_jet_pt : sub_pt) {
      for (auto& sub_const_pt_init : const_sub_pt_init) {
        for (auto& sub_const_pt_match : const_sub_pt_match) {
          for (auto& sub_const_eta : const_eta) {
              
            // get a few needed parameters
            double R = sub_fj_def.R();
            double bkg_R = bkg_definition.R();
            double jet_eta_max = sub_const_eta - R;
            double bkg_jet_eta_max = sub_const_eta - bkg_R;
              
            // build constituent & jet selectors
            fastjet::Selector init_const_selector = fastjet::SelectorAbsRapMax(sub_const_eta)
                                                  && fastjet::SelectorPtMin(sub_const_pt_init);
            fastjet::Selector match_const_selector = fastjet::SelectorAbsRapMax(sub_const_eta)
                                                   && fastjet::SelectorPtMin(sub_const_pt_match);
            fastjet::Selector init_jet_selector = fastjet::SelectorAbsRapMax(jet_eta_max)
                                                && fastjet::SelectorPtMin(sub_jet_pt);
            fastjet::Selector match_jet_selector = fastjet::SelectorIdentity();
            fastjet::Selector bkg_selector = fastjet::SelectorAbsRapMax(sub_const_eta - bkg_R)
                                          && (!fastjet::SelectorNHardest(2));
              
            fastjet::GhostedAreaSpec ghost_def(jet_eta_max + 2 * R, ghost_repeat, ghost_area,
                                               grid_scatter, pt_scatter, mean_ghost_pt);
            fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts,
                                              ghost_def);
            
            fastjet::GhostedAreaSpec bkg_ghost_def(bkg_jet_eta_max + 2 * bkg_R, ghost_repeat, ghost_area,
                                                   grid_scatter, pt_scatter, mean_ghost_pt);
            fastjet::AreaDefinition bkg_area_def(fastjet::active_area_explicit_ghosts,
                                                 bkg_ghost_def);
              
            // create the initial JetDef & the
            // matched JetDef def
            JetDef init_def(sub_fj_def, area_def, bkg_definition, bkg_area_def, bkg_selector);
            init_def.SetConstituentSelector(init_const_selector);
            init_def.SetJetSelector(init_jet_selector);
            JetDef match_def(sub_fj_def, area_def, bkg_definition, bkg_area_def, bkg_selector);
            match_def.SetConstituentSelector(match_const_selector);
            match_def.SetJetSelector(match_jet_selector);
            
            // create the MatchDef
            sub_matchdefs.insert(std::make_unique<MatchDef>(init_def, match_def));
            
          }
        }
      }
    }
  }
  
  // create the dijet definitions
  for (auto& lead : lead_matchdefs) {
    for (auto& sub : sub_matchdefs) {
      
      // make sure we're looking at similar jets
      if (lead->InitialJetDef().jet_algorithm() != sub->InitialJetDef().jet_algorithm())
        continue;
      // make sure that the pt lead > pt sub
      if (SelectorPtMinLessThan(lead->InitialJetDef().JetSelector(), sub->InitialJetDef().JetSelector()))
        continue;
      // if forcce_constituent_pt_equality is on, force
      // pt to be equal. Same for eta
      if (force_constituent_pt_equality) {
        double lead_pt_init = ExtractDoubleFromSelector(lead->InitialJetDef().ConstituentSelector(), "pt >=");
        double sub_pt_init = ExtractDoubleFromSelector(sub->InitialJetDef().ConstituentSelector(), "pt >=");
        if (lead_pt_init != sub_pt_init)
          continue;
        double lead_pt_match = ExtractDoubleFromSelector(lead->MatchedJetDef().ConstituentSelector(), "pt >=");
        double sub_pt_match = ExtractDoubleFromSelector(sub->MatchedJetDef().ConstituentSelector(), "pt >=");
        if (lead_pt_match != sub_pt_match)
          continue;
      }
      if (force_constituent_eta_equality) {
        double lead_eta = ExtractDoubleFromSelector(lead->InitialJetDef().ConstituentSelector(), "|rap| <=");
        double sub_eta = ExtractDoubleFromSelector(sub->InitialJetDef().ConstituentSelector(), "|rap| <=");
        if (lead_eta != sub_eta)
          continue;
      }
      
      auto tmp = std::make_unique<DijetDefinition>(lead.get(), sub.get(), 0.4);
      
      std::string key = MakeKeyFromDijetDefinition(*tmp);
      dijet_defs.insert({key, std::move(tmp)});
      keys.insert(key);
    }
  }
  
}

void DijetMatrix::CheckToUpdate() {
  if (dijet_defs.size() != 0) {
    ClearDijetDefs();
  }
}

void DijetMatrix::InitializeEmptyFields() {
  if (const_eta.empty())
    const_eta.insert(1.0);
  if (const_lead_pt_init.empty())
    const_lead_pt_init.insert(2.0);
  if (const_lead_pt_match.empty())
    const_lead_pt_match.insert(0.2);
  if (const_sub_pt_init.empty())
    const_sub_pt_init.insert(2.0);
  if (const_sub_pt_match.empty())
    const_sub_pt_match.insert(0.2);
  if (jet_algorithm.empty())
    jet_algorithm.insert(fastjet::antikt_algorithm);
  if (lead_pt.empty())
    lead_pt.insert(20.0);
  if (lead_R.empty())
    lead_R.insert(0.4);
  if (sub_pt.empty())
    sub_pt.insert(10.0);
  if (sub_R.empty())
    sub_R.insert(0.4);
}

std::vector<fastjet::JetDefinition> DijetMatrix::FillLeadJetDefinitions() {
  std::vector<fastjet::JetDefinition> ret;
  for (auto& alg : jet_algorithm)
    for (auto& R : lead_R)
      ret.push_back(fastjet::JetDefinition(alg, R, scheme, strategy));
  return ret;
}

std::vector<fastjet::JetDefinition> DijetMatrix::FillSubJetDefinitions() {
  std::vector<fastjet::JetDefinition> ret;
  for (auto& alg : jet_algorithm)
    for (auto& R : sub_R)
      ret.push_back(fastjet::JetDefinition(alg, R, scheme, strategy));
  return ret;
}
