// match_def.cc

#include "jet_analysis/dijet_worker/dijet_matrix/match_def.hh"

bool MatchDef::EquivalentCluster(const MatchDef& rhs) const {
  return rhs.InitialJetDef().EquivalentCluster(InitialJetDef()) &&
         rhs.MatchedJetDef().EquivalentCluster(InitialJetDef());
}
