// dijet_matrix.hh

// the DijetMatrix defines a set of parameters
// that defines a dijet definition (as used by
// the dijet imbalance measurement). The user
// can specify a set of values for any/each parameter,
// and the DijetMatrix will create a DijetDefinition
// for each valid set of the parameters.

// .... obviously, having multiple values in, say, 10
// parameters is a bad idea. (even if only two values are
// specified for each parameter, the total number of
// DijetDefinitions would be massive).

#ifndef DIJET_MATRIX_HH
#define DIJET_MATRIX_HH

#include <memory>
#include <set>
#include <unordered_map>
#include <vector>

#include "jet_analysis/dijet_worker/dijet_definition.hh"
#include "jet_analysis/dijet_worker/match_def.hh"
#include "jet_analysis/dijet_worker/jet_def.hh"

#include "fastjet/JetDefinition.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/Selector.hh"

class DijetMatrix {
  
public:
  
  // default constructor populates all fields that aren't filled
  // via InitializeEmptyFields()
  DijetMatrix();
  
  // single entry construction
  DijetMatrix(fastjet::JetAlgorithm jet_alg_in,
              double lead_jet_pt_in = 20.0,
              double lead_jet_R_in = 0.4,
              double sub_jet_pt_in = 10.0,
              double sub_jet_R_in = 0.4,
              double const_lead_pt_init_in = 2.0,
              double const_lead_pt_match_in = 0.2,
              double const_sub_pt_init_in = 2.0,
              double const_sub_pt_match_in = 0.2,
              double eta_in = 1.0);
  
  
  // set construction
  DijetMatrix(std::set<fastjet::JetAlgorithm> jet_alg_in,
              std::set<double> lead_pt_in,
              std::set<double> lead_R_in,
              std::set<double> sub_pt_in,
              std::set<double> sub_R_in,
              std::set<double> const_lead_pt_init_in,
              std::set<double> const_lead_pt_match_in,
              std::set<double> const_sub_pt_init_in,
              std::set<double> const_sub_pt_match_in,
              std::set<double> eta_in);
  
  DijetMatrix(const DijetMatrix& rhs);
  
  virtual ~DijetMatrix() { };
  
  // removes all parameters
  void Clear();
  // removes only the DijetDefinitions, if Initialize()
  // has been called
  void ClearDijetDefs();
  
  // used to initialize the matrix.
  void Initialize();
  
  // get the number of dijet definitions currently stored
  // (will be zero if no initialization has occured)
  inline std::size_t Size() {return dijet_defs.size();}
  
  // get access to the dijet definitions
  std::unordered_map<std::string, std::shared_ptr<DijetDefinition>> DijetDefinitions() const
                                                                           {return dijet_defs;}
  
  // get the keys to the map
  const std::set<std::string>& Keys() const {return keys;}
  
  // The Matrix will build all logically consistent di-jet pairs
  // (i.e. if pt sub > pt lead, it will be ignored, that sort of
  // thing). Adding extra parameters will increase the effective
  // run time considerably, especially as the set of non-unique
  // parameters increases.
  void AddConstituentEta(double eta);
  void AddConstituentEta(std::set<double> eta);
  
  void AddConstituentLeadInitialPt(double pt);
  void AddConstituentLeadInitialPt(std::set<double> pt);
  
  void AddConstituentLeadMatchPt(double pt);
  void AddConstituentLeadMatchPt(std::set<double> pt);
  
  void AddConstituentSubInitialPt(double pt);
  void AddConstituentSubInitialPt(std::set<double> pt);
  
  void AddConstituentSubMatchPt(double pt);
  void AddConstituentSubMatchPt(std::set<double> pt);
  
  void AddJetAlgorithm(fastjet::JetAlgorithm alg);
  void AddJetAlgorithm(std::set<fastjet::JetAlgorithm> alg);
  
  void AddLeadJetPt(double pt);
  void AddLeadJetPt(std::set<double> pt);
  
  void AddLeadJetR(double R);
  void AddLeadJetR(std::set<double> R);
  
  void AddSubJetPt(double pt);
  void AddSubJetPt(std::set<double> pt);
  
  void AddSubJetR(double R);
  void AddSubJetR(std::set<double> R);
  
  // options to change the default fastjet settings
  // for area/bkg estimation
  inline void SetClusterStrategy(fastjet::Strategy strat)  {strategy = strat;}
  inline void SetRecombinationScheme(fastjet::RecombinationScheme schm)
                                                           {scheme = schm;}
  inline void SetAreaType(fastjet::AreaType type)          {area_type = type;}
  inline void SetGhostRepeat(int repeat)                   {ghost_repeat = repeat;}
  inline void SetGhostArea(double area)                    {ghost_area = area;}
  inline void SetGridScatter(double scatter)               {grid_scatter = scatter;}
  inline void SetPtScatter(double scatter)                 {pt_scatter = scatter;}
  inline void SetMeanGhostPt(double mean)                  {mean_ghost_pt = mean;}
  inline void SetBkgDefinition(fastjet::JetDefinition def) {bkg_definition = def;}
  
  // access to the internally stored sets
  inline const std::set<double>& ConstituentEta() const              {return const_eta;}
  inline const std::set<double>& LeadConstituentInitPt() const       {return const_lead_pt_init;}
  inline const std::set<double>& LeadConstituentMatchPt() const      {return const_lead_pt_match;}
  inline const std::set<double>& SubConstituentInitPt() const        {return const_sub_pt_init;}
  inline const std::set<double>& SubConstituentMatchPt() const       {return const_sub_pt_match;}
  inline const std::set<fastjet::JetAlgorithm>& JetAlgorithm() const {return jet_algorithm;}
  inline const std::set<double>& LeadJetPt() const                   {return lead_pt;}
  inline const std::set<double>& LeadJetR() const                    {return lead_R;}
  inline const std::set<double>& SubJetPt() const                    {return sub_pt;}
  inline const std::set<double>& SubJetR() const                     {return sub_R;}
  
  inline fastjet::Strategy ClusterStrategy() const                {return strategy;}
  inline fastjet::RecombinationScheme RecombinationScheme() const {return scheme;}
  inline fastjet::AreaType AreaType() const                       {return area_type;}
  
  
protected:
  
  // used internally to update the dijet definitions, if there has
  // been a parameter updated after initialization
  void CheckToUpdate();
  
  // used internally to make sure there can be at least one valid
  // dijet definition when Initialize is called. Otherwise, it sets
  // default values for missing fields
  void InitializeEmptyFields();
  
  // used during initialization
  std::vector<fastjet::JetDefinition> FillLeadJetDefinitions();
  std::vector<fastjet::JetDefinition> FillSubJetDefinitions();
  
  std::set<double> const_eta;
  std::set<double> const_lead_pt_init;
  std::set<double> const_lead_pt_match;
  std::set<double> const_sub_pt_init;
  std::set<double> const_sub_pt_match;
  
  // allow different radii for leading & subleading
  // jets, but do not allow different jet algorithms
  std::set<fastjet::JetAlgorithm> jet_algorithm;
  std::set<double> lead_pt;
  std::set<double> lead_R;
  std::set<double> sub_pt;
  std::set<double> sub_R;
  
  // details of fastjet that one might want to change from default,
  // for one reason or another
  
  // recombination strategy, probably shouldn't change. Default is set
  // to 'best' which tells fastjet to choose the most efficient
  fastjet::Strategy strategy;
  
  // recombination scheme
  fastjet::RecombinationScheme scheme;
  
  // the type of area estimation done, default is active area with
  // explicit ghosts
  fastjet::AreaType area_type;
  
  // ghost repeating, default = 1
  int ghost_repeat;
  
  // ghost area (smaller = more ghosts per unit area)
  double ghost_area;
  
  // for completeness, these are included - can look up in fastjet
  // manual for details.
  double grid_scatter;
  double pt_scatter;
  double mean_ghost_pt;
  
  // background jet definition
  fastjet::JetDefinition bkg_definition;
  
  // populated by the DijetMatrix
  std::unordered_map<std::string, std::shared_ptr<DijetDefinition>> dijet_defs;
  std::set<std::shared_ptr<MatchDef>> lead_matchdefs;
  std::set<std::shared_ptr<MatchDef>> sub_matchdefs;
  
  // keys for the dictionary of DijetDefintions are identifier
  // strings describing the detailed settings of each definition.
  // the keys are the minimum length needed to be unique to each
  // dijet definition
  std::set<std::string> keys;
};

#endif // DIJET_MATRIX_HH
