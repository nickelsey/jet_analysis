// dijet_definition.cc

#include "jet_analysis/dijet_worker/dijet_matrix/dijet_definition.hh"

bool DijetDefinition::EquivalentCluster(const DijetDefinition& rhs) const {
  return lead->EquivalentCluster(*rhs.lead) && sub->EquivalentCluster(*rhs.sub);
}
