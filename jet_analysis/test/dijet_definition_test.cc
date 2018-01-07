#include "jet_analysis/dijet_worker/jet_def.hh"
#include "jet_analysis/dijet_worker/match_def.hh"
#include "jet_analysis/dijet_worker/dijet_definition.hh"

int main() {

  DijetDefinition default_def;

  if (*default_def.lead != MatchDef() ||
      *default_def.sub != MatchDef())
    return 1;
  
  if (default_def.IsValid() != false ||
      default_def.DoMatching() != false)
    return 1;
  
  JetDef valid_jetdef(fastjet::antikt_algorithm, 0.4);
  std::shared_ptr<MatchDef> valid_matchdef =
      std::make_shared<MatchDef>(valid_jetdef, valid_jetdef);
  DijetDefinition valid_dijetdef(valid_matchdef, valid_matchdef);
  
  if (*valid_dijetdef.lead != *valid_matchdef ||
      *valid_dijetdef.sub != *valid_matchdef)
    return 1;
  
  if (valid_dijetdef.IsValid() != true ||
      valid_dijetdef.DoMatching() != true)
    return 1;
  
  JetDef default_jetdef;
  std::shared_ptr<MatchDef> valid_matchdef_no_matching =
      std::make_shared<MatchDef>(valid_jetdef, default_jetdef);
  DijetDefinition valid_dijetdef_no_matching(valid_matchdef_no_matching, valid_matchdef_no_matching);
  if (valid_dijetdef_no_matching.IsValid() != true ||
      valid_dijetdef_no_matching.DoMatching() != false)
    return 1;
  
  DijetDefinition tmp(valid_dijetdef_no_matching);
  if (tmp != valid_dijetdef_no_matching)
    return 1;
  
  return 0;
}
