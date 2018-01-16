// pythia_aj_efficiency.cc

#include "jet_analysis/util/arg_helper.hh"
#include "jet_analysis/dijet_worker/dijet_worker.hh"
#include "jet_analysis/util/pt_distribution.hh"
#include "jet_analysis/util/string_util.hh"
#include "jet_analysis/util/vector_conversion.hh"

#include <string>
#include <iostream>
#include <random>

#include "Pythia8/Pythia.h"

#include "fastjet/PseudoJet.hh"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TStyle.h"
#include "TFile.h"

#include "boost/filesystem.hpp"

using std::string;
struct Options {
  string pt_low  = "25";    /* lower pt */
  string pt_high = "100";   /* upper pt */
  string out_dir = "tmp";   /* location to place output */
  string name    = "job";   /* name of job */
  int id         = 0;       /* which job, if being submitted concurrently */
  int ev         = 2e3;     /* number of events per efficiency bin */
  int n_bkg      = 400;     /* number of background particles per event */
};

int main(int argc, char* argv[]) {
  
  constexpr double pi = 3.14159265359;
  
  // Histograms will calculate gaussian errors
  // -----------------------------------------
  TH1::SetDefaultSumw2( );
  TH2::SetDefaultSumw2( );
  TH3::SetDefaultSumw2( );
  
  // set drawing preferences for histograms and graphs
  gStyle->SetOptStat(false);
  gStyle->SetOptFit(false);
  gStyle->SetOptTitle(1);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetHatchesSpacing(1.0);
  gStyle->SetHatchesLineWidth(2);
  
  Options opts;
  for (int i = 1; i < argc; ++i) {
    if (ParseStrFlag(string(argv[i]), "--name", &opts.name) ||
        ParseStrFlag(string(argv[i]), "--outDir", &opts.out_dir) ||
        ParseStrFlag(string(argv[i]), "--ptLow", &opts.pt_low) ||
        ParseStrFlag(string(argv[i]), "--ptHigh", &opts.pt_high) ||
        ParseIntFlag(string(argv[i]), "--eventsPerBin", &opts.ev) ||
        ParseIntFlag(string(argv[i]), "--nbkg", &opts.n_bkg) ||
        ParseIntFlag(string(argv[i]), "--id", &opts.id)) continue;
    std::cerr << "unknown command line option: " << argv[i] << std::endl;
    return 1;
  }
  
  // build output directory if it doesn't exist, using boost::filesystem
  if (opts.out_dir.empty())
    opts.out_dir = "tmp";
  boost::filesystem::path dir(opts.out_dir.c_str());
  boost::filesystem::create_directories(dir);
  
  TFile out(MakeString(opts.out_dir, "/", opts.name,
                       opts.id, ".root").c_str(), "RECREATE");
  
  // number of events per efficiency bin
  size_t events_per_bin = opts.ev;
  
  // create the generator
  // --------------------
  string pt_low_string = "PhaseSpace:pTHatMin = 20.0";
  if (!opts.pt_low.empty())
    pt_low_string = "PhaseSpace:pTHatMin = " + opts.pt_low;
  string pt_high_string = "PhaseSpace:pTHatMax = 100.0";
  if (!opts.pt_high.empty())
    pt_high_string = "PhaseSpace:pTHatMax = " + opts.pt_high;
  string qcd_string = "HardQCD:all = on";
  string e_cm_string = "Beams:eCM = 200";
  Pythia8::Pythia pythia;
  pythia.readString(qcd_string);
  pythia.readString(e_cm_string);
  pythia.readString(pt_low_string);
  pythia.readString(pt_high_string);
  pythia.readString("Random:setSeed = On");
  pythia.readString("Random:seed = 0");
  
  pythia.init();
  
  // create generators for a random background
  std::uniform_real_distribution<double> prob_dis(0, 1.0);
  std::uniform_real_distribution<double> eta_dis(-1.0, 1.0);
  std::uniform_real_distribution<double> phi_dis(0.0, 2*pi);
  pt_distribution<double> pt_dis;
  
  // create the generator
  // standard mersenne_twister_engine
  std::random_device device;
  std::mt19937 gen(device());
  
  // dijet worker
  DijetWorker worker;
  worker.Initialize();
  auto keys = worker.Keys();
  
  // Use star kinematics
  fastjet::Selector particle_selector = fastjet::SelectorPtMin(0.2) && fastjet::SelectorAbsRapMax(1.0);
  
  // define the efficiencies we'll use
  std::vector<double> efficiencies{1.0, 0.9, 0.8, 0.7, 0.6, 0.5};
  std::vector<string> eff_string{"100%", "90%", "80%", "70%", "60%", "50%"};
  
  // define histograms for the efficiencies
  std::vector<TH1D*> h_Aj;
  std::vector<TH1D*> h_nPart_before;
  std::vector<TH1D*> h_nPart_after;
  
  for (size_t eff = 0; eff < efficiencies.size(); ++eff) {
    
    h_Aj.push_back(new TH1D(MakeString("aj_", eff).c_str(), "aj", 30, 0, 0.9));
    h_nPart_after.push_back(new TH1D(MakeString("npa_", eff).c_str(), "npart after", 100, 0.5, 600.5));
    h_nPart_before.push_back(new TH1D(MakeString("npb_", eff).c_str(), "npart before", 100, 0.5, 600.5));
    
    double efficiency = efficiencies[eff];
    
    // generate events_per_bin events per efficiency bin
    for (size_t ev = 0; ev < events_per_bin * (eff + 1); ++ev) {
      
      if (!pythia.next())
        continue;
      
      std::vector<fastjet::PseudoJet> particles(opts.n_bkg);
      
      // add a random background of opts.n_bkg particles
      for (size_t i = 0; i < opts.n_bkg; ++i) {
        fastjet::PseudoJet tmp;
        tmp.reset_PtYPhiM(pt_dis(gen), eta_dis(gen), phi_dis(gen), 0);
        tmp.set_user_index(10);
        particles[i] = tmp;
      }
      
      particles = particle_selector(particles);
      
      ConvertPythiaToPseudoJet(pythia, particles);
      
      std::vector<fastjet::PseudoJet> input(particles.size());
      
      size_t accepted = 0;
      for (auto vec : particles) {
        if (prob_dis(gen) >= efficiency)
          continue;
        input[accepted] = vec;
        ++accepted;
      }
      input.resize(accepted);
      
      auto worker_out = worker.Run(input);
      
      for (auto def_out : worker_out) {
        if (!def_out.second.found_match)
          continue;
        double Aj = (def_out.second.lead_hard.pt() - def_out.second.sublead_hard.pt()) /
                    (def_out.second.lead_hard.pt() + def_out.second.sublead_hard.pt());
        h_Aj[eff]->Fill(Aj);
        h_nPart_before[eff]->Fill(particles.size());
        h_nPart_after[eff]->Fill(input.size());
      }
    }
  }
  
  out.Write();
  out.Close();
  
  return 0;
}
