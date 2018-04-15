// run7_eff.hh


#ifndef RUN7_EFF_HH
#define RUN7_EFF_HH

#include <vector>
#include <string>

#include "TH2.h"
#include "TF2.h"
#include "TFile.h"

class Run7Eff {
public:
  Run7Eff();
  Run7Eff(std::string filename);
  
  ~Run7Eff();
  
  void loadFile(std::string filename);
  
  double AuAuEff(double pt, double eta, int cent);
  double AuAuEff020Avg(double pt, double eta);
  double ppEff(double pt, double eta);
  
private:
  
  void loadCurves();
  TF2* GetEffY06();
  TF2* GetEffY04(int cb);
  
  TFile* file;
  
  //centrality bins: 0 is 0-5%, 1 is 5-10%, 2 is 10-20%
  TF2* effY04[3]; // Run 4 parameterization
  TH2D* effY07pteta[3]; // Run 7 HT / Run 4 MB pt-eta map, used for pt <= 1.5 GeV/c
  TH1D* effY07eta[3]; // Run 7 HT / Run 4 MB eta map, used for pt > 1.5 GeV/c
  TF2* effY06; // Run 6 parameterization
  
  double maxPt;
  double maxPtpp;
  
};

#endif // RUN7_EFF_HH
