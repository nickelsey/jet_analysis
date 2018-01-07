// dijet_definition.hh

// A DijetDefinition contains a MatchDef for
// both the leading & subleading jet

// direct access is given to the MatchDefs, since
// this is just a container

#ifndef DIJET_DEFINITION_HH
#define DIJET_DEFINITION_HH

#include <memory>

#include "jet_analysis/dijet_worker/match_def.hh"

class DijetDefinition {
public:
  DijetDefinition() : lead(std::make_shared<MatchDef>()), sub(std::make_shared<MatchDef>()), dPhi(0.4) { };
  DijetDefinition(const DijetDefinition& rhs) :
                  lead(rhs.lead), sub(rhs.sub), dPhi(rhs.dPhi) { };
  DijetDefinition(std::shared_ptr<MatchDef> lead, std::shared_ptr<MatchDef> sub, double phi = 0.4) :
                  lead(lead), sub(sub), dPhi(phi) { };
  
  virtual ~DijetDefinition() { };
  
  inline bool IsValid()    const {return lead && sub && lead->IsValid() && sub->IsValid();}
  inline bool DoMatching() const {return lead && sub && lead->CanMatch() && sub->CanMatch();}
  
  std::shared_ptr<MatchDef> lead;
  std::shared_ptr<MatchDef> sub;
  
  double dPhi;
};

// tests for logical inequality, not pointer
// inequality - so if the pointers differ, but
// the underlying MatchDefs are equivalent, will
// return true
inline bool operator==(const DijetDefinition& lhs, const DijetDefinition& rhs) {
  if (lhs.dPhi != rhs.dPhi)
    return false;
  if (lhs.lead == nullptr && lhs.sub == nullptr &&
      rhs.lead == nullptr && rhs.sub == nullptr)
    return true;
  if ((lhs.lead == nullptr && rhs.lead != nullptr) ||
      (lhs.lead != nullptr && rhs.lead == nullptr) ||
      (lhs.sub == nullptr  && rhs.sub != nullptr) ||
      (lhs.sub != nullptr  && rhs.sub == nullptr))
    return false;
  if (lhs.lead != nullptr &&
      *lhs.lead != *rhs.lead)
    return false;
  if (lhs.sub != nullptr &&
      *lhs.sub != *rhs.sub)
    return false;
  return true;
}
inline bool operator!=(const DijetDefinition& lhs, const DijetDefinition& rhs) {
  return !(lhs == rhs);
}

#endif // DIJET_DEFINITION_HH
