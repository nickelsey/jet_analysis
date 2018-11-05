// run11_eff.hh


#ifndef RUN4_EFF_HH
#define RUN4_EFF_HH

#include <vector>
#include <string>

#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"

class Run11Eff {
public:
  Run11Eff(std::string filename);
  
  ~Run11Eff() {};

	bool LoadFile(std::string filename);
  
  double AuAuEff(double pt, double eta, int cent);
  double AuAuEff(double pt, int cent);
  
private:

	TFile* input_;
  std::vector<TH2D*> curves_;
  std::vector<TH1D*> curves_1d_;
  
  double max_pt_;
  
};


#endif // RUN11_EFF_HH