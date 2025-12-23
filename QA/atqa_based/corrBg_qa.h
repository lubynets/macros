//
// Created by oleksii on 23.12.2025.
//
#ifndef ANALYSISTREEQA_CORRBG_QA_H
#define ANALYSISTREEQA_CORRBG_QA_H

#include <cstdint>
#include <map>
#include <string>
#include <vector>
#include <utility>

/// Channels taken from here: https://github.com/AliceO2Group/O2Physics/blob/87be5da87be8bcef56dccf64b5d960e7f1b7545d/PWGHF/Core/DecayChannels.h#L61-L95
/// @brief 3-prong candidates: main channels
enum DecayChannelMain : int8_t {
  // D+
  DplusToPiKPi = 1,    // π+ K− π+
  DplusToPiKPiPi0 = 2, // π+ K− π+ π0
  DplusToPiPiPi = 3,   // π+ π− π+
  DplusToPiKK = 4,     // π+ K− K+
  // Ds+
  DsToPiKK = 5,      // π+ K− K+
  DsToPiKKPi0 = 6,   // π+ K− K+ π0
  DsToPiPiK = 7,     // π+ π− K+
  DsToPiPiPi = 8,    // π+ π− π+
  DsToPiPiPiPi0 = 9, // π+ π− π+ π0
  // D*+
  DstarToPiKPi = 10,       // π+ K− π+ (from [(D0 → π+ K−) π+])
  DstarToPiKPiPi0 = 11,    // π+ K− π+ π0
  DstarToPiKPiPi0Pi0 = 12, // π+ K− π+ π0 π0
  DstarToPiKK = 13,        // π+ K− K+
  DstarToPiKKPi0 = 14,     // π+ K− K+ π0
  DstarToPiPiPi = 15,      // π+ π− π+
  DstarToPiPiPiPi0 = 16,   // π+ π− π+ π0
  // Λc+
  LcToPKPi = 17,    // p K− π+
  LcToPKPiPi0 = 18, // p K− π+ π0
  LcToPPiPi = 19,   // p π− π+
  LcToPKK = 20,     // p K− K+
  // Ξc+
  XicToPKPi = 21,  // p K− π+
  XicToPKK = 22,   // p K− K+
  XicToSPiPi = 23, // Σ+ π− π+
  //
  NChannelsMain = XicToSPiPi // last channel
};

std::map<int, float> BRsDplus_PYTHIA = {
  std::make_pair(DplusToPiKPi,    0.28),
  std::make_pair(DplusToPiKPiPi0, 0.05),
  std::make_pair(DplusToPiPiPi,   0.34),
  std::make_pair(DplusToPiKK,     0.33)
};

std::map<int, float> BRsDs_PYTHIA = {
  std::make_pair(DsToPiKK,        0.20),
  std::make_pair(DsToPiKKPi0,     0.05),
  std::make_pair(DsToPiPiK,       0.25),
  std::make_pair(DsToPiPiPi,      0.25),
  std::make_pair(DsToPiPiPiPi0,   0.25)
};

std::map<int, float> BRsD0_PYTHIA = {
  std::make_pair(DstarToPiKPi,       /* 0.6770000 * FOR D*: TO BE ADDED IN THE WEIGHT CALCULATION  */ 0.20),
  std::make_pair(DstarToPiKPiPi0,    /* 0.6770000 * FOR D*: TO BE ADDED IN THE WEIGHT CALCULATION  */ 0.20),
  std::make_pair(DstarToPiKPiPi0Pi0, /* 0.6770000 * FOR D*: TO BE ADDED IN THE WEIGHT CALCULATION  */ 0.0000000000001 /*not simulated*/),
  std::make_pair(DstarToPiKK,        /* 0.6770000 * FOR D*: TO BE ADDED IN THE WEIGHT CALCULATION  */ 0.20),
  std::make_pair(DstarToPiKKPi0,     /* 0.6770000 * FOR D*: TO BE ADDED IN THE WEIGHT CALCULATION  */ 0.0000000000001 /*not simulated*/),
  std::make_pair(DstarToPiPiPi,      /* 0.6770000 * FOR D*: TO BE ADDED IN THE WEIGHT CALCULATION  */ 0.20),
  std::make_pair(DstarToPiPiPiPi0,   /* 0.6770000 * FOR D*: TO BE ADDED IN THE WEIGHT CALCULATION  */ 0.20)
};

std::map<int, float> BRsLc_PYTHIA = {
  std::make_pair(LcToPKPi,    0.20),
  std::make_pair(LcToPKPiPi0, 0.05),
  std::make_pair(LcToPPiPi,   0.25),
  std::make_pair(LcToPKK,     0.25)
};

std::map<int, float> BRsXic_PYTHIA = {
  std::make_pair(XicToPKPi,   0.34),
  std::make_pair(XicToPKK,    0.33),
  std::make_pair(XicToSPiPi,  0.33)
};

/// @brief D+ BR values from PDG (https://pdg.lbl.gov/2025/tables/contents_tables_mesons.html)
std::map<int, float> BRsDplus_PDG = {
  std::make_pair(DplusToPiKPi,    0.0752 /*(K pi)_s_wave pi*/),
  std::make_pair(DplusToPiKPiPi0, 0.0625),
  std::make_pair(DplusToPiPiPi,   0.00084 /*rho0 pi*/ + 0.000454 /*f2(1270) pi*/ + 0.00201 /*pi (pi pi)_s_wave*/),
  std::make_pair(DplusToPiKK,     0.00249 /*K*0 K*/ + 0.00182 /*K*0(1430) K*/ + 0.00269 /*phi pi*/ + 0.00268 /*non-resonant, put as phi pi*/)
};
/// @brief Ds+ BR values from PDG (https://pdg.lbl.gov/2025/tables/contents_tables_mesons.html)
std::map<int, float> BRsDs_PDG = {
  std::make_pair(DsToPiKK,        0.0225 /*phi pi*/ + 0.02261 /*K* K*/),
  std::make_pair(DsToPiKKPi0,     0.0559 /*phi rho*/),
  std::make_pair(DsToPiPiK,       0.00168 /*K* pi*/ + 0.0012 /*f0(1370) K*/ + 0.00218 /*rho0 K*/ + 0.00099 /*non-resonant*/),
  std::make_pair(DsToPiPiPi,      0.00014 /*rho0 pi*/ + 0.00142 /*f2(1270) pi*/ + 0.00923 /*pi(pi pi)_s_wave*/),
  std::make_pair(DsToPiPiPiPi0,   0.01686 /*eta pi+*/)
};
/// @brief D0 BR values from PDG (https://pdg.lbl.gov/2025/tables/contents_tables_mesons.html) --> BR(D* -> D0pi) factorizes when computing the scaling factor ==> it can be ignored
std::map<int, float> BRsD0_PDG = {
  std::make_pair(DstarToPiKPi,       /* 0.6770000 * FOR D*: TO BE ADDED IN THE WEIGHT CALCULATION  */ 0.03945),
  std::make_pair(DstarToPiKPiPi0,    /* 0.6770000 * FOR D*: TO BE ADDED IN THE WEIGHT CALCULATION  */ 0.015 /*non-resonant*/ + 0.112 /*rho+ K-*/ + 0.0195 /*K*0 pi0*/ + 0.0231 /*K*- pi+*/),
  std::make_pair(DstarToPiKPiPi0Pi0, /* 0.6770000 * FOR D*: TO BE ADDED IN THE WEIGHT CALCULATION  */ 0.0886),
  std::make_pair(DstarToPiKK,        /* 0.6770000 * FOR D*: TO BE ADDED IN THE WEIGHT CALCULATION  */ 0.00408),
  std::make_pair(DstarToPiKKPi0,     /* 0.6770000 * FOR D*: TO BE ADDED IN THE WEIGHT CALCULATION  */ 0.00342 /*total*/),
  std::make_pair(DstarToPiPiPi,      /* 0.6770000 * FOR D*: TO BE ADDED IN THE WEIGHT CALCULATION  */ 0.001453),
  std::make_pair(DstarToPiPiPiPi0,   /* 0.6770000 * FOR D*: TO BE ADDED IN THE WEIGHT CALCULATION  */ 0.0101 /*rho+ pi-*/ + 0.00013 /*non-resonant*/)
};
/// @brief Lc+ BR values from PDG (https://pdg.lbl.gov/2025/tables/contents_tables_baryons.html)
std::map<int, float> BRsLc_PDG = {
  std::make_pair(LcToPKPi,    0.0635),
  std::make_pair(LcToPKPiPi0, 0.046 /*non-resonant*/),
  std::make_pair(LcToPPiPi,   0.00467 /*total*/),
  std::make_pair(LcToPKK,     0.00105 /*p phi*/ * 0.499 /*phi->K+K-*/)
};
/// @brief Xic+ BR values values from PDG (https://pdg.lbl.gov/2025/tables/contents_tables_baryons.html)
std::map<int, float> BRsXic_PDG = {
  std::make_pair(XicToPKPi,   0.0062 /*total*/),
  std::make_pair(XicToPKK,    0.00012 /*p phi*/ * 0.499 /*phi->K+K-*/),
  std::make_pair(XicToSPiPi,  0.014)
};

/// Scaling factor in PYTHIA for the D+ decays, from the fact that PYTHIA normalizes all BRs to sum at 1
const float pythiaBRScalingFactorDPlus = 1./( BRsDplus_PYTHIA[DplusToPiKPi] + BRsDplus_PYTHIA[DplusToPiKPiPi0] + BRsDplus_PYTHIA[DplusToPiPiPi] + BRsDplus_PYTHIA[DplusToPiKK] );
/// Scaling factor in PYTHIA for the Ds+ decays, from the fact that PYTHIA normalizes all BRs to sum at 1
const float pythiaBRScalingFactorDs = 1./( BRsDs_PYTHIA[DsToPiKK] + BRsDs_PYTHIA[DsToPiKKPi0] + BRsDs_PYTHIA[DsToPiPiK] + BRsDs_PYTHIA[DsToPiPiPi] + BRsDs_PYTHIA[DsToPiPiPiPi0]);
/// Scaling factor in PYTHIA for the D*+ decays, from the fact that PYTHIA normalizes all BRs to sum at 1
const float pythiaBRScalingFactorDstar = 1./( BRsD0_PYTHIA[DstarToPiKPi] + BRsD0_PYTHIA[DstarToPiKPiPi0] + BRsD0_PYTHIA[DstarToPiKPiPi0Pi0] + BRsD0_PYTHIA[DstarToPiKK] + BRsD0_PYTHIA[DstarToPiKKPi0] + BRsD0_PYTHIA[DstarToPiPiPi] + BRsD0_PYTHIA[DstarToPiPiPiPi0]);
/// Scaling factor in PYTHIA for the Lc+ decays, from the fact that PYTHIA normalizes all BRs to sum at 1
const float pythiaBRScalingFactorLc = 1./( BRsLc_PYTHIA[LcToPKPi] + BRsLc_PYTHIA[LcToPKPiPi0] + BRsLc_PYTHIA[LcToPPiPi] + BRsLc_PYTHIA[LcToPKK] );
/// Scaling factor in PYTHIA for the Xic+ decays, from the fact that PYTHIA normalizes all BRs to sum at 1
const float pythiaBRScalingFactorXic = 1./( BRsXic_PYTHIA[XicToPKPi] + BRsXic_PYTHIA[XicToPKK] + BRsXic_PYTHIA[XicToSPiPi] );

struct Mother {
  int id_;
  std::string name_;
  std::string greek_name_;
  float pythia_br_scaling_factor_;
};

enum MotherParticle : int {
  Dplus = 0,
  Ds,
  Dstar,
  Lc,
  Xic,
  NMotherParticles
};

const std::vector<Mother> Mothers {
  {Dplus, "Dplus", "D^{+}",           pythiaBRScalingFactorDPlus},
  {Ds,    "Ds",    "D_{s}^{+}",       pythiaBRScalingFactorDs   },
  {Dstar, "Dstar", "D^{*+}",          pythiaBRScalingFactorDstar},
  {Lc,    "Lc",    "#Lambda_{c}^{+}", pythiaBRScalingFactorLc   },
  {Xic,   "Xic",   "#Xi_{c}^{+}",     pythiaBRScalingFactorXic  }
};

struct Decay {
  int id_;
  Mother mother_;
  std::string daughters_;
  float br_pythia_;
  float br_pdg_;
  std::string greek_daughters_;
  int fill_style_;
  bool is_bg_;
};

const std::vector<Decay> Decays {
  {DplusToPiKPi,       Mothers[Dplus], "PiKPi",       BRsDplus_PYTHIA[DplusToPiKPi],    BRsDplus_PDG[DplusToPiKPi],    "#pi^{+}K^{#minus}#pi^{+}",               3004, true },
  {DplusToPiKPiPi0,    Mothers[Dplus], "PiKPiPi0",    BRsDplus_PYTHIA[DplusToPiKPiPi0], BRsDplus_PDG[DplusToPiKPiPi0], "#pi^{+}K^{#minus}#pi^{+}#pi^{0}",        3005, true },
  {DplusToPiPiPi,      Mothers[Dplus], "PiPiPi",      BRsDplus_PYTHIA[DplusToPiPiPi],   BRsDplus_PDG[DplusToPiPiPi],   "#pi^{+}#pi^{#minus}#pi^{+}",             3006, true },
  {DplusToPiKK,        Mothers[Dplus], "PiKK",        BRsDplus_PYTHIA[DplusToPiKK],     BRsDplus_PDG[DplusToPiKK],     "K^{+}K^{#minus}#pi^{+}",                 3007, true },

  {DsToPiKK,           Mothers[Ds],    "PiKK",        BRsDs_PYTHIA[DsToPiKK],           BRsDs_PDG[DsToPiKK],           "K^{+}K^{#minus}#pi^{+}",                 3004, true },
  {DsToPiKKPi0,        Mothers[Ds],    "PiKKPi0",     BRsDs_PYTHIA[DsToPiKKPi0],        BRsDs_PDG[DsToPiKKPi0],        "K^{+}K^{#minus}#pi^{+}#pi^{0}",          3005, true },
  {DsToPiPiK,          Mothers[Ds],    "PiPiK",       BRsDs_PYTHIA[DsToPiPiK],          BRsDs_PDG[DsToPiPiK],          "#pi^{+}#pi^{#minus}K^{+}",               3006, true },
  {DsToPiPiPi,         Mothers[Ds],    "PiPiPi",      BRsDs_PYTHIA[DsToPiPiPi],         BRsDs_PDG[DsToPiPiPi],         "#pi^{+}#pi^{#minus}#pi^{+}",             3007, true },
  {DsToPiPiPiPi0,      Mothers[Ds],    "PiPiPiPi0",   BRsDs_PYTHIA[DsToPiPiPiPi0],      BRsDs_PDG[DsToPiPiPiPi0],      "#pi^{+}#pi^{#minus}#pi^{+}#pi^{0}",      3008, true },

  {DstarToPiKPi,       Mothers[Dstar], "PiKPi",       BRsD0_PYTHIA[DstarToPiKPi],       BRsD0_PDG[DstarToPiKPi],       "#pi^{+}K^{#minus}#pi^{+}",               3004, true },
  {DstarToPiKPiPi0,    Mothers[Dstar], "PiKPiPi0",    BRsD0_PYTHIA[DstarToPiKPiPi0],    BRsD0_PDG[DstarToPiKPiPi0],    "#pi^{+}K^{#minus}#pi^{+}#pi^{0}",        3005, true },
  {DstarToPiKPiPi0Pi0, Mothers[Dstar], "PiKPiPi0Pi0", BRsD0_PYTHIA[DstarToPiKPiPi0Pi0], BRsD0_PDG[DstarToPiKPiPi0Pi0], "#pi^{+}K^{#minus}#pi^{+}#pi^{0}#pi^{0}", 3006, true },
  {DstarToPiKK,        Mothers[Dstar], "PiKK",        BRsD0_PYTHIA[DstarToPiKK],        BRsD0_PDG[DstarToPiKK],        "K^{+}K^{#minus}#pi^{+}",                 3007, true },
  {DstarToPiKKPi0,     Mothers[Dstar], "PiKKPi0",     BRsD0_PYTHIA[DstarToPiKKPi0],     BRsD0_PDG[DstarToPiKKPi0],     "K^{+}K^{#minus}#pi^{+}#pi^{0}",          3008, true },
  {DstarToPiPiPi,      Mothers[Dstar], "PiPiPi",      BRsD0_PYTHIA[DstarToPiPiPi],      BRsD0_PDG[DstarToPiPiPi],      "#pi^{+}#pi^{#minus}#pi^{+}",             3009, true },
  {DstarToPiPiPiPi0,   Mothers[Dstar], "PiPiPiPi0",   BRsD0_PYTHIA[DstarToPiPiPiPi0],   BRsD0_PDG[DstarToPiPiPiPi0],   "#pi^{+}#pi^{#minus}#pi^{+}#pi^{0}",      3010, true },

  {LcToPKPi,           Mothers[Lc],    "PKPi",        BRsLc_PYTHIA[LcToPKPi],           BRsLc_PDG[LcToPKPi],           "pK^{#minus}#pi^{+}",                     3001, true },
  {LcToPKPi,           Mothers[Lc],    "PKPi",        BRsLc_PYTHIA[LcToPKPi],           BRsLc_PDG[LcToPKPi],           "pK^{#minus}#pi^{+}",                     3002, false},
  {LcToPKPiPi0,        Mothers[Lc],    "PKPiPi0",     BRsLc_PYTHIA[LcToPKPiPi0],        BRsLc_PDG[LcToPKPiPi0],        "pK^{#minus}#pi^{+}#pi^{0}",              3004, true },
  {LcToPPiPi,          Mothers[Lc],    "PPiPi",       BRsLc_PYTHIA[LcToPPiPi],          BRsLc_PDG[LcToPPiPi],          "p#pi^{#minus}#pi^{+}",                   3005, true },
  {LcToPKK,            Mothers[Lc],    "PKK",         BRsLc_PYTHIA[LcToPKK],            BRsLc_PDG[LcToPKK],            "pK^{#minus}K^{+}",                       3006, true },

  {XicToPKPi,          Mothers[Xic],   "PKPi",        BRsXic_PYTHIA[XicToPKPi],         BRsXic_PDG[XicToPKPi],         "pK^{#minus}#pi^{+}",                     3004, true },
  {XicToPKK,           Mothers[Xic],   "PKK",         BRsXic_PYTHIA[XicToPKK],          BRsXic_PDG[XicToPKK],          "pK^{#minus}K^{+}",                       3005, true },
  {XicToSPiPi,         Mothers[Xic],   "SPiPi",       BRsXic_PYTHIA[XicToSPiPi],        BRsXic_PDG[XicToSPiPi],        "#Sigma^{+}#pi^{#minus}#pi^{+}",          3006, true }
};

#endif//ANALYSISTREEQA_CORRBG_QA_H
