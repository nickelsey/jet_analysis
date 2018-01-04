// vector_conversion
// functions to easily convert between
// TLorentzVector, fastjet::PseudoJet, etc

#ifndef VECTOR_CONVERSION_HH
#define VECTOR_CONVERSION_HH

#include "fastjet/PseudoJet.hh"
#include "TStarJetVectorContainer.h"
#include "TStarJetVector.h"

#include <random>

// converts a TStarJetVectorContainer of TStarJetVectors
// into an std::vector of fastjet::PseudoJets
void ConvertTStarJetVector(TStarJetVectorContainer<TStarJetVector>* container, std::vector<fastjet::PseudoJet>& particles) {
  
  // Transform TStarJetVectors into (FastJet) PseudoJets
  // ---------------------------------------------------
  TStarJetVector* sv;
  for(int i = 0; i < container->GetEntries() ; ++i) {
    sv = container->Get(i);
    
    fastjet::PseudoJet tmpPJ = fastjet::PseudoJet(*sv);
    tmpPJ.set_user_index(sv->GetCharge());
    particles.push_back(tmpPJ);
  }
}

#endif // VECTOR_CONVERSION_HH
