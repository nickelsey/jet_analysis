// run14_eff.cc

#include <iostream>

#include "run14_eff.hh"

Run14Eff::Run14Eff() : file(nullptr), maxPt(4.5) {
  
}

Run14Eff::Run14Eff(std::string filename) : file(nullptr), maxPt(4.5) {
  loadFile(filename);
}

Run14Eff::~Run14Eff() {
  
}

void Run14Eff::loadFile(std::string filename) {
  if (file)
    file->Close();
  file = new TFile(filename.c_str(), "READ");
  loadCurves();
}

double Run14Eff::AuAuEff(double pt, double eta, int cent, double zdcrate) {
  if (pt > maxPt)
    pt = maxPt;
  
  int zdcBin;
  if (zdcrate/1000.0 < 33.0)
    zdcBin = 0;
  else if (zdcrate/1000.0 < 66.0)
    zdcBin = 1;
  else
    zdcBin = 2;
  
  std::cout << "getting pt: " << pt << std::endl;
  std::cout << "zdc bin: " << zdcBin << std::endl;
  std::cout << "cent bin: " << cent << std::endl;
  std::cout << "max zdc size: " << curves.size() << std::endl;
  std::cout << "max cent size: " << curves[0].size() << std::endl;
  
  int bin = curves[zdcBin][cent]->FindBin(pt, eta);
  return curves[zdcBin][cent]->GetBinContent(bin);
  
}

void Run14Eff::loadCurves(int nBinsZDC, int nBinsCent) {
  curves.clear();
  for (int i = 0; i < nBinsZDC; ++i) {
    curves.push_back(std::vector<TH2D*>());
    for (int j = 0; j < nBinsCent; ++j) {
      std::string name = "efficiency_lumi_" + std::to_string(i)
                       + "_cent_" + std::to_string(j);
      curves[i].push_back((TH2D*) file->Get(name.c_str()));
    }
  }
}
