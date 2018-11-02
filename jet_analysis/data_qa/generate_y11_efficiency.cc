// generate_y11_efficiency.cc

#include "jet_analysis/util/common.hh"
#include "jet_analysis/util/string_util.hh"

#include <string>
#include <iostream>
#include <vector>

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

DEFINE_string(name, "y11_efficiency", "output file name");
DEFINE_string(input, "run11_embed.root", "input Run 11 embedding file");
DEFINE_string(outdir, "y11_efficiency", "output directory");

int main(int argc, char* argv[]) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);
  
  // build output directory if it doesn't exist, using boost::filesystem
  if (FLAGS_outdir.empty()) FLAGS_outdir = "tmp";
  boost::filesystem::path dir(FLAGS_outdir.c_str());
  boost::filesystem::create_directories(dir);

	TFile input(FLAGS_input.c_str(), "READ");
	TH3D* mc = (TH3D*) input.Get("mc");
	TH3D* match = (TH3D*) input.Get("match");

	// create output file from the given directory, name & id
  string outfile_name = FLAGS_outdir + "/" + FLAGS_name + ".root";
  TFile out(outfile_name.c_str(), "RECREATE");

	std::vector<TH2D*> mc_projection;
	std::vector<TH2D*> match_projection;

	for (int i = 0; i < 9; ++i) {
		mc->GetZaxis()->SetRange(i+1, i+1);
		match->GetZaxis()->SetRange(i+1, i+1);

		std::string mc_name = "mc" + std::to_string(i);
		std::string match_name = "cent" + std::to_string(i);

		mc_projection.push_back((TH2D*) mc->Project3D("yx"));
		mc_projection[i]->SetName(mc_name.c_str());
		match_projection.push_back((TH2D*) match->Project3D("yx"));
		match_projection[i]->SetName(match_name.c_str());
		mc_projection[i]->RebinX(5);
		match_projection[i]->RebinX(5);

		match_projection[i]->Divide(mc_projection[i]);
	}

	for (auto& h : match_projection)
		h->Write();
	out.Close();

	return 0;
}