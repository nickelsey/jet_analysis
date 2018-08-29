// jet_def.cc

#include "jet_analysis/dijet_worker/dijet_matrix/jet_def.hh"

// JetDef class imp

JetDef::JetDef() :
  fastjet::JetDefinition(),
  const_selector(fastjet::SelectorIdentity()),
  jet_selector(fastjet::SelectorIdentity()),
  area_def(fastjet::AreaDefinition(fastjet::invalid_area, fastjet::GhostedAreaSpec())),
  bkg_jet_def(fastjet::JetDefinition()),
  bkg_area_def(fastjet::AreaDefinition(fastjet::invalid_area, fastjet::GhostedAreaSpec())),
  bkg_selector(!fastjet::SelectorNHardest(2)) { }

JetDef::JetDef(fastjet::JetAlgorithm alg, double R,
               fastjet::AreaDefinition area_def,
               fastjet::JetDefinition bkg_jet_def,
               fastjet::AreaDefinition bkg_area_def,
               fastjet::Selector bkg_selector) :
  fastjet::JetDefinition(alg, R),
  const_selector(fastjet::SelectorIdentity()),
  jet_selector(fastjet::SelectorIdentity()),
  area_def(area_def),
  bkg_jet_def(bkg_jet_def),
  bkg_area_def(bkg_area_def),
  bkg_selector(bkg_selector) { }

JetDef::JetDef(fastjet::JetAlgorithm alg, double R,
               double xtra_param,
               fastjet::AreaDefinition area_def,
               fastjet::JetDefinition bkg_jet_def,
               fastjet::AreaDefinition bkg_area_def,
               fastjet::Selector bkg_selector) :
  fastjet::JetDefinition(alg, R, xtra_param),
  const_selector(fastjet::SelectorIdentity()),
  jet_selector(fastjet::SelectorIdentity()),
  area_def(area_def),
  bkg_jet_def(bkg_jet_def),
bkg_area_def(bkg_area_def),
  bkg_selector(bkg_selector) { }

JetDef::JetDef(fastjet::JetDefinition jet_def_in,
               fastjet::AreaDefinition area_def,
               fastjet::JetDefinition bkg_jet_def,
               fastjet::AreaDefinition bkg_area_def,
               fastjet::Selector bkg_selector) :
  fastjet::JetDefinition(jet_def_in),
  const_selector(fastjet::SelectorIdentity()),
  jet_selector(fastjet::SelectorIdentity()),
  area_def(area_def),
  bkg_jet_def(bkg_jet_def),
  bkg_area_def(bkg_area_def),
  bkg_selector(bkg_selector) { }

JetDef::JetDef(const JetDef& rhs) :
  fastjet::JetDefinition(rhs.jet_algorithm(), rhs.R(),
                         rhs.recombination_scheme(),
                         rhs.strategy()),
  const_selector(rhs.const_selector),
  jet_selector(rhs.jet_selector),
  area_def(rhs.area_def),
  bkg_jet_def(rhs.bkg_jet_def),
  bkg_area_def(rhs.bkg_area_def),
  bkg_selector(rhs.bkg_selector) { }

bool JetDef::EquivalentCluster(const JetDef& rhs) const {
  // essentially == operator, without jet selector comparison
  if (jet_algorithm() != rhs.jet_algorithm() ||
      R() != rhs.R() ||
      recombination_scheme() != rhs.recombination_scheme() ||
      strategy() != rhs.strategy() ||
      ConstituentSelector().description() != rhs.ConstituentSelector().description() ||
      AreaDefinition().area_type() != rhs.AreaDefinition().area_type() ||
      AreaDefinition().ghost_spec().ghost_maxrap() != rhs.AreaDefinition().ghost_spec().ghost_maxrap() ||
      AreaDefinition().ghost_spec().ghost_area() != rhs.AreaDefinition().ghost_spec().ghost_area() ||
      AreaDefinition().ghost_spec().grid_scatter() != rhs.AreaDefinition().ghost_spec().grid_scatter() ||
      AreaDefinition().ghost_spec().pt_scatter() != rhs.AreaDefinition().ghost_spec().pt_scatter() ||
      AreaDefinition().ghost_spec().mean_ghost_pt() != rhs.AreaDefinition().ghost_spec().mean_ghost_pt() ||
      AreaDefinition().ghost_spec().repeat() != rhs.AreaDefinition().ghost_spec().repeat() ||
      BackgroundAreaDef().area_type() != rhs.BackgroundAreaDef().area_type() ||
      BackgroundAreaDef().ghost_spec().ghost_maxrap() != rhs.BackgroundAreaDef().ghost_spec().ghost_maxrap() ||
      BackgroundAreaDef().ghost_spec().ghost_area() != rhs.BackgroundAreaDef().ghost_spec().ghost_area() ||
      BackgroundAreaDef().ghost_spec().grid_scatter() != rhs.BackgroundAreaDef().ghost_spec().grid_scatter() ||
      BackgroundAreaDef().ghost_spec().pt_scatter() != rhs.BackgroundAreaDef().ghost_spec().pt_scatter() ||
      BackgroundAreaDef().ghost_spec().mean_ghost_pt() != rhs.BackgroundAreaDef().ghost_spec().mean_ghost_pt() ||
      BackgroundAreaDef().ghost_spec().repeat() != rhs.BackgroundAreaDef().ghost_spec().repeat() ||
      BackgroundJetDef().jet_algorithm() != rhs.BackgroundJetDef().jet_algorithm() ||
      BackgroundJetDef().R() != rhs.BackgroundJetDef().R() ||
      BackgroundJetDef().recombination_scheme() != rhs.BackgroundJetDef().recombination_scheme() ||
      BackgroundJetDef().strategy() != rhs.BackgroundJetDef().strategy() ||
      BackgroundSelector().description() != rhs.BackgroundSelector().description())
    return false;
  return true;
}



