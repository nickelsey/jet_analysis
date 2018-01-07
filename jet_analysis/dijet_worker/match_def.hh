// match_def.hh

#ifndef MATCH_DEF_HH
#define MATCH_DEF_HH

#include "jet_analysis/dijet_worker/jet_def.hh"

class MatchDef {
public:
  
  // By default, both initial & match definitions
  // are in invalid states to prevent misunderstanding
  MatchDef() : initial(), matched(), R(0.0) { };
  
  // if not specified, dR is set to the minimum of the
  // two JetDef radii
  MatchDef(const JetDef& init, const JetDef& match) :
    initial(init), matched(match) {
    R = initial.R() < matched.R() ? initial.R() : matched.R();
  }
  
  // can specify a non-default radius
  MatchDef(const JetDef& init, const JetDef& match,
           double dR) : initial(init), matched(match), R(dR) { };
  
  virtual ~MatchDef() { }
  
  inline JetDef InitialJetDef() const {return initial;}
  inline JetDef MatchedJetDef() const {return matched;}
  
  inline double  dR() const {return R;}
  
  void SetInitialJetDef(const JetDef& def) {initial = def;}
  void SetMatchedJetDef(const JetDef& def) {matched = def;}
  
  void SetdR(double dR_in) {R = dR_in;}
  
  // matching is done when a both JetDefs are in valid states
  inline bool IsValid()  const  {return initial.IsValid();}
  inline bool CanMatch() const  {return IsValid() && matched.IsValid();}
  
private:
  
  JetDef initial;
  JetDef matched;
  
  double R;
  
};

inline bool operator==(const MatchDef& lhs, const MatchDef& rhs) {
  if (lhs.InitialJetDef() != rhs.InitialJetDef() ||
      lhs.MatchedJetDef() != rhs.MatchedJetDef() ||
      lhs.dR() != rhs.dR())
    return false;
  return true;
}

inline bool operator!=(const MatchDef& lhs, const MatchDef& rhs) {
  return !(lhs == rhs);
}

#endif // MATCH_DEF_HH
