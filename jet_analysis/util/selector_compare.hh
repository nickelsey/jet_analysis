// selector_compare.hh

// used to compare selectors based on their string description
// to keep things as simple as possible, I have decided that
// when applicable, these will act like operators
// (eg, SelectorPtMinGreaterThan(lhs, rhs) will return true if lhs' pt cut
// is greater than rhs's).
// To avoid throwing exceptions, it will be assumed that if a selector
// does not contain a specific cut, then its value is taken to be zero


#ifndef SELECTOR_COMPARE_HH
#define SELECTOR_COMPARE_HH

#include <sstream>

#include "fastjet/Selector.hh"

double ExtractDoubleFromSelector(const fastjet::Selector& sel, std::string descriptor) {
  
  std::string sel_str = sel.description();
  std::size_t length = descriptor.size();
  std::size_t pos    = sel_str.find(descriptor);
  
  if ( pos == std::string::npos)
    return 0.0;
  
  std::istringstream stream(std::string(sel_str.begin() + pos + length, sel_str.end()));
  double ret;
  stream >> std::skipws >> ret;
  return ret;
}

bool SelectorPtMinGreaterThan(const fastjet::Selector& lhs, const fastjet::Selector& rhs) {
  std::string descriptor = "pt >=";
  
  double lhs_pt = ExtractDoubleFromSelector(lhs, descriptor);
  double rhs_pt = ExtractDoubleFromSelector(rhs, descriptor);

  return lhs_pt > rhs_pt;
}

bool SelectorPtMinGreaterEqualThan(const fastjet::Selector& lhs, const fastjet::Selector& rhs) {
  std::string descriptor = "pt >=";
  
  double lhs_pt = ExtractDoubleFromSelector(lhs, descriptor);
  double rhs_pt = ExtractDoubleFromSelector(rhs, descriptor);
  
  return lhs_pt >= rhs_pt;
}

bool SelectorPtMinLessThan(const fastjet::Selector& lhs, const fastjet::Selector& rhs) {
  std::string descriptor = "pt >=";
  
  double lhs_pt = ExtractDoubleFromSelector(lhs, descriptor);
  double rhs_pt = ExtractDoubleFromSelector(rhs, descriptor);
  
  return lhs_pt < rhs_pt;
}

bool SelectorPtMinLessEqualThan(const fastjet::Selector& lhs, const fastjet::Selector& rhs) {
  std::string descriptor = "pt >=";
  
  double lhs_pt = ExtractDoubleFromSelector(lhs, descriptor);
  double rhs_pt = ExtractDoubleFromSelector(rhs, descriptor);
  
  return lhs_pt <= rhs_pt;
}

bool SelectorPtMaxGreaterThan(const fastjet::Selector& lhs, const fastjet::Selector& rhs) {
  std::string descriptor = "pt <=";
  
  double lhs_pt = ExtractDoubleFromSelector(lhs, descriptor);
  double rhs_pt = ExtractDoubleFromSelector(rhs, descriptor);
  
  return lhs_pt > rhs_pt;
}

bool SelectorPtMaxGreaterEqualThan(const fastjet::Selector& lhs, const fastjet::Selector& rhs) {
  std::string descriptor = "pt <=";
  
  double lhs_pt = ExtractDoubleFromSelector(lhs, descriptor);
  double rhs_pt = ExtractDoubleFromSelector(rhs, descriptor);
  
  return lhs_pt >= rhs_pt;
}

bool SelectorPtMaxLessThan(const fastjet::Selector& lhs, const fastjet::Selector& rhs) {
  std::string descriptor = "pt <=";
  
  double lhs_pt = ExtractDoubleFromSelector(lhs, descriptor);
  double rhs_pt = ExtractDoubleFromSelector(rhs, descriptor);
  
  return lhs_pt < rhs_pt;
}

bool SelectorPtMaxLessEqualThan(const fastjet::Selector& lhs, const fastjet::Selector& rhs) {
  std::string descriptor = "pt <=";
  
  double lhs_pt = ExtractDoubleFromSelector(lhs, descriptor);
  double rhs_pt = ExtractDoubleFromSelector(rhs, descriptor);
  
  return lhs_pt <= rhs_pt;
}


#endif //SELECTOR_COMPARE_HH
