// run14_eff.hh


#ifndef RUN14_EFF_HH
#define RUN14_EFF_HH

#include <vector>
#include <string>
#include <memory>

#include "TH2.h"
#include "TF2.h"
#include "TFile.h"

class Run14Eff {
public:
  Run14Eff(std::string filename = "submit/y14_effic.root");
  
  ~Run14Eff();
  
  void loadFile(std::string filename);
  
  double AuAuEff(double pt, double eta, int cent, double zdcrate);
  double pp6Eff(double pt, double eta);
  double ratio(double pt, double eta, int cent, double zdcrate);
  
private:
  
  void loadCurves(int nBinsZDC = 3, int nBinsCent = 11);
  
  TFile* file;
  
  std::vector<std::vector<TH2D*>> curves;
  std::shared_ptr<TF2> effY06; // Run 6 parameterization
                
  double maxPt;
  
};


#endif // RUN14_EFF_HH
