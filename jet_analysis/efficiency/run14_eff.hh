// run14_eff.hh


#ifndef RUN14_EFF_HH
#define RUN14_EFF_HH

#include <vector>
#include <string>
#include <memory>

#include "TH2.h"
#include "TF2.h"
#include "TFile.h"

enum class TrackingUnc {
  NONE = 0,
  POSITIVE = 1,
  NEGATIVE = -1
};

class Run14Eff {
public:
  
  Run14Eff(std::string filename = "submit/y14_effic_dca3.root");
  
  ~Run14Eff();
  
  void loadFile(std::string filename, int nBinsZDC = 3, int nBinsCent = 16);
  
  double AuAuEff(double pt, double eta, int cent, double zdcrate);
  double pp6Eff(double pt, double eta);
  double ratio(double pt, double eta, int cent, double zdcrate);
  double ratioUncertainty(double pt, double eta, int cent, double zdcrate);
  
  void setSystematicUncertainty(TrackingUnc sys = TrackingUnc::NONE) {sys_ = sys;}
  TrackingUnc SystematicUncertainty() const {return sys_;}
  
  void setAuAuUncertainty(double unc) {auau_u_ = unc;}
  double AuAuUncertainty() const {return auau_u_;}
  void setPPUncertainty(double unc) {pp_u_ = unc;}
  double PPUncertainty() const {return pp_u_;}
  void setCentUncertainty(double unc) {cent_u_ = unc;}
  double CentUncertainty() const {return cent_u_;}
  
private:
  
  void loadCurves(int nBinsZDC = 3, int nBinsCent = 16);
  
  std::shared_ptr<TFile> file;
  
  std::vector<std::vector<TH2D*>> curves;
  std::shared_ptr<TF2> effY06; // Run 6 parameterization
                
  double max_pt_;
  
  TrackingUnc sys_;
  double auau_u_;
  double pp_u_;
  double cent_u_;
  
};


#endif // RUN14_EFF_HH
