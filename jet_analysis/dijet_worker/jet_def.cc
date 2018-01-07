// jet_def.cc

#include "jet_analysis/dijet_worker/jet_def.hh"

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



