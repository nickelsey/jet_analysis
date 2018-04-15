// run14_eff.hh


#ifndef RUN14_EFF_HH
#define RUN14_EFF_HH

#include <vector>
#include <string>

#include "TH2.h"
#include "TFile.h"

class Run14Eff {
public:
  Run14Eff();
  Run14Eff(std::string filename);
  
  ~Run14Eff();
  
  void loadFile(std::string filename);
  
  double AuAuEff(double pt, double eta, int cent, double zdcrate);
  
private:
  
  void loadCurves(int nBinsZDC = 3, int nBinsCent = 11);
  
  TFile* file;
  
  std::vector<std::vector<TH2D*>> curves;
  
  double maxPt;
  
};


#endif // RUN14_EFF_HH
