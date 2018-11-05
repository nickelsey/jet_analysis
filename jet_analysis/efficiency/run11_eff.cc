#include "jet_analysis/efficiency/run11_eff.hh"
#include <iostream>
Run11Eff::Run11Eff(std::string file) {
	LoadFile(file);
}

bool Run11Eff::LoadFile(std::string file) {
	input_ = new TFile(file.c_str(), "READ");
	for (int i = 0; i < 9; ++i) {
		std::string curvename = "cent" + std::to_string(i);
		curves_.push_back((TH2D*) input_->Get(curvename.c_str()));
	}
	return true;
}

double Run11Eff::AuAuEff(double pt, double eta, int cent) {
	if (cent < 0 || cent > 8)
		return 1.0;

	if (pt < 0.2)
		return 1.0;

	if (pt > max_pt_)
		pt = max_pt_;
	
	if (fabs(eta) > 1.0)
		return 1.0;
	std::cout << "finding bins" << std::endl;
	int binx = curves_.at(cent)->GetXaxis()->FindBin(pt);
	int biny = curves_.at(cent)->GetYaxis()->FindBin(eta);
	std::cout << "finding content" << std::endl;
	double efficiency = curves_.at(cent)->GetBinContent(binx, biny);
	return efficiency;
}