#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TColor.h"

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <map>

template <typename I>
std::string to_string_fixed_digits(I value,int digitsCount);
//__________________________________________________________
//__________________________________________________________
//__________________________________________________________
std::vector<std::string> GetDFNames(const std::string& fileName);
//__________________________________________________________
//__________________________________________________________
//__________________________________________________________
void setRowCols(const int size, int& rows, int& cols, int& w, int& h);
//__________________________________________________________
//__________________________________________________________
//__________________________________________________________
/// Channels taken from here: https://github.com/AliceO2Group/O2Physics/blob/master/PWGHF/Core/DecayChannels.h#L61-L95 (8th July 2025)
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
/*[x]*/  DstarToPiKPi = 10,       // π+ K− π+ (from [(D0 → π+ K−) π+])
/*[x]*/  DstarToPiKPiPi0 = 11,    // π+ K− π+ π0
/*[ ]*/  DstarToPiKPiPi0Pi0 = 12, // π+ K− π+ π0 π0
/*[x]*/  DstarToPiKK = 13,        // π+ K− K+
/*[ ]*/  DstarToPiKKPi0 = 14,     // π+ K− K+ π0
/*[x]*/  DstarToPiPiPi = 15,      // π+ π− π+
/*[x]*/  DstarToPiPiPiPi0 = 16,   // π+ π− π+ π0
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

/// @brief D+ BR values used in MC (https://github.com/AliceO2Group/O2DPG/blob/master/MC/config/PWGHF/pythia8/generator/pythia8_charmhadronic_with_decays_Mode2_CorrBkg.cfg)
std::map<int, float> BRsDplus_PYTHIA = {
    std::make_pair(DplusToPiKPi,    0.28),
    std::make_pair(DplusToPiKPiPi0, 0.05),
    std::make_pair(DplusToPiPiPi,   0.34),
    std::make_pair(DplusToPiKK,     0.33)
};
/// @brief Ds+ BR values used in MC (https://github.com/AliceO2Group/O2DPG/blob/master/MC/config/PWGHF/pythia8/generator/pythia8_charmhadronic_with_decays_Mode2_CorrBkg.cfg)
std::map<int, float> BRsDs_PYTHIA = {
    std::make_pair(DsToPiKK,        0.20),
    std::make_pair(DsToPiKKPi0,     0.05),
    std::make_pair(DsToPiPiK,       0.25),
    std::make_pair(DsToPiPiPi,      0.25),
    std::make_pair(DsToPiPiPiPi0,   0.25)
};
/// @brief D0 BR values used in MC (https://github.com/AliceO2Group/O2DPG/blob/master/MC/config/PWGHF/pythia8/generator/pythia8_charmhadronic_with_decays_Mode2_CorrBkg.cfg)
std::map<int, float> BRsD0_PYTHIA = {
    std::make_pair(DstarToPiKPi,       /* 0.6770000 * FOR D*: TO BE ADDED IN THE WEIGHT CALCULATION  */ 0.20),
    std::make_pair(DstarToPiKPiPi0,    /* 0.6770000 * FOR D*: TO BE ADDED IN THE WEIGHT CALCULATION  */ 0.20),
    std::make_pair(DstarToPiKPiPi0Pi0, /* 0.6770000 * FOR D*: TO BE ADDED IN THE WEIGHT CALCULATION  */ 0.0000000000001 /*not simulated*/),
    std::make_pair(DstarToPiKK,        /* 0.6770000 * FOR D*: TO BE ADDED IN THE WEIGHT CALCULATION  */ 0.20),
    std::make_pair(DstarToPiKKPi0,     /* 0.6770000 * FOR D*: TO BE ADDED IN THE WEIGHT CALCULATION  */ 0.0000000000001 /*not simulated*/),
    std::make_pair(DstarToPiPiPi,      /* 0.6770000 * FOR D*: TO BE ADDED IN THE WEIGHT CALCULATION  */ 0.20),
    std::make_pair(DstarToPiPiPiPi0,   /* 0.6770000 * FOR D*: TO BE ADDED IN THE WEIGHT CALCULATION  */ 0.20)
};
/// @brief Lc+ BR values used in MC (https://github.com/AliceO2Group/O2DPG/blob/master/MC/config/PWGHF/pythia8/generator/pythia8_charmhadronic_with_decays_Mode2_CorrBkg.cfg)
std::map<int, float> BRsLc_PYTHIA = {
    std::make_pair(LcToPKPi,    0.20),
    std::make_pair(LcToPKPiPi0, 0.05),
    std::make_pair(LcToPPiPi,   0.25),
    std::make_pair(LcToPKK,     0.25)
};
/// @brief Xic+ BR values used in MC (https://github.com/AliceO2Group/O2DPG/blob/master/MC/config/PWGHF/pythia8/generator/pythia8_charmhadronic_with_decays_Mode2_CorrBkg.cfg)
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
float pythiaBRScalingFactorDPlus = 1./( BRsDplus_PYTHIA[DplusToPiKPi] + BRsDplus_PYTHIA[DplusToPiKPiPi0] + BRsDplus_PYTHIA[DplusToPiPiPi] + BRsDplus_PYTHIA[DplusToPiKK] );
/// Scaling factor in PYTHIA for the Ds+ decays, from the fact that PYTHIA normalizes all BRs to sum at 1
float pythiaBRScalingFactorDs = 1./( BRsDs_PYTHIA[DsToPiKK] + BRsDs_PYTHIA[DsToPiKKPi0] + BRsDs_PYTHIA[DsToPiPiK] + BRsDs_PYTHIA[DsToPiPiPi] + BRsDs_PYTHIA[DsToPiPiPiPi0]);
/// Scaling factor in PYTHIA for the D*+ decays, from the fact that PYTHIA normalizes all BRs to sum at 1
float pythiaBRScalingFactorDstar = 1./( BRsD0_PYTHIA[DstarToPiKPi] + BRsD0_PYTHIA[DstarToPiKPiPi0] + BRsD0_PYTHIA[DstarToPiKPiPi0Pi0] + BRsD0_PYTHIA[DstarToPiKK] + BRsD0_PYTHIA[DstarToPiKKPi0] + BRsD0_PYTHIA[DstarToPiPiPi] + BRsD0_PYTHIA[DstarToPiPiPiPi0]);
/// Scaling factor in PYTHIA for the Lc+ decays, from the fact that PYTHIA normalizes all BRs to sum at 1
float pythiaBRScalingFactorLc = 1./( BRsLc_PYTHIA[LcToPKPi] + BRsLc_PYTHIA[LcToPKPiPi0] + BRsLc_PYTHIA[LcToPPiPi] + BRsLc_PYTHIA[LcToPKK] );
/// Scaling factor in PYTHIA for the Xic+ decays, from the fact that PYTHIA normalizes all BRs to sum at 1
float pythiaBRScalingFactorXic = 1./( BRsXic_PYTHIA[XicToPKPi] + BRsXic_PYTHIA[XicToPKK] + BRsXic_PYTHIA[XicToSPiPi] );

struct Mother {
  int id_;
  std::string name_;
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

std::vector<Mother> Mothers {
  {Dplus, "Dplus", pythiaBRScalingFactorDPlus},
  {Ds,    "Ds",    pythiaBRScalingFactorDs},
  {Dstar, "Dstar", pythiaBRScalingFactorDstar},
  {Lc,    "Lc",    pythiaBRScalingFactorLc},
  {Xic,   "Xic",   pythiaBRScalingFactorXic}
};

struct Decay {
  int id_;
  Mother mother_;
  std::string daughters_;
  float br_pythia_;
  float br_pdg_;
  std::string formula_;
  int fill_style_;
};
std::vector<Decay> Decays {
  {DplusToPiKPi,        Mothers[Dplus],    "PiKPi",        BRsDplus_PYTHIA[DplusToPiKPi],        BRsDplus_PDG[DplusToPiKPi],        "D^{+}#rightarrowK^{#minus}#pi^{+}#pi^{+}",                3004},
  {DplusToPiKPiPi0,     Mothers[Dplus],    "PiKPiPi0",     BRsDplus_PYTHIA[DplusToPiKPiPi0],     BRsDplus_PDG[DplusToPiKPiPi0],     "D^{+}#rightarrowK^{#minus}#pi^{+}#pi^{+}#pi^{0}",         3005},
  {DplusToPiPiPi,       Mothers[Dplus],    "PiPiPi",       BRsDplus_PYTHIA[DplusToPiPiPi],       BRsDplus_PDG[DplusToPiPiPi],       "D^{+}#rightarrow#pi^{#minus}#pi^{+}#pi^{+}",              3006},
  {DplusToPiKK,         Mothers[Dplus],    "PiKK",         BRsDplus_PYTHIA[DplusToPiKK],         BRsDplus_PDG[DplusToPiKK],         "D^{+}#rightarrow#pi^{+}K^{#minus}K^{+}",                  3007},

  {DsToPiKK,            Mothers[Ds],       "PiKK",         BRsDs_PYTHIA[DsToPiKK],               BRsDs_PDG[DsToPiKK],               "D_{s}^{+}#rightarrowK^{#minus}#pi^{+}#pi^{+}",             3004},
  {DsToPiKKPi0,         Mothers[Ds],       "PiKKPi0",      BRsDs_PYTHIA[DsToPiKKPi0],            BRsDs_PDG[DsToPiKKPi0],            "D_{s}^{+}#rightarrowK^{#minus}#pi^{+}#pi^{+}#pi^{0}",      3005},
  {DsToPiPiK,           Mothers[Ds],       "PiPiK",        BRsDs_PYTHIA[DsToPiPiK],              BRsDs_PDG[DsToPiPiK],              "D_{s}^{+}#rightarrowK^{+}#pi^{+}#pi^{#minus}",             3006},
  {DsToPiPiPi,          Mothers[Ds],       "PiPiPi",       BRsDs_PYTHIA[DsToPiPiPi],             BRsDs_PDG[DsToPiPiPi],             "D_{s}^{+}#rightarrow#pi^{+}#pi^{#minus}#pi^{+}",           3007},
  {DsToPiPiPiPi0,       Mothers[Ds],       "PiPiPiPi0",    BRsDs_PYTHIA[DsToPiPiPiPi0],          BRsDs_PDG[DsToPiPiPiPi0],          "D_{s}^{+}#rightarrow#pi^{+}#pi^{#minus}#pi^{+}#pi^{0}",    3008},

  {DstarToPiKPi,        Mothers[Dstar],    "PiKPi",        BRsD0_PYTHIA[DstarToPiKPi],            BRsD0_PDG[DstarToPiKPi],          "D^{*+}#rightarrow#pi^{+}K^{#minus}#pi^{+}",                3004},
  {DstarToPiKPiPi0,     Mothers[Dstar],    "PiKPiPi0",     BRsD0_PYTHIA[DstarToPiKPiPi0],         BRsD0_PDG[DstarToPiKPiPi0],       "D^{*+}#rightarrow#pi^{+}K^{#minus}#pi^{+}#pi^{0}",         3005},
  {DstarToPiKPiPi0Pi0,  Mothers[Dstar],    "PiKPiPi0Pi0",  BRsD0_PYTHIA[DstarToPiKPiPi0Pi0],      BRsD0_PDG[DstarToPiKPiPi0Pi0],    "D^{*+}#rightarrow#pi^{+}K^{#minus}#pi^{+}#pi^{0}#pi^{0}",  3006},
  {DstarToPiKK,         Mothers[Dstar],    "PiKK",         BRsD0_PYTHIA[DstarToPiKK],             BRsD0_PDG[DstarToPiKK],           "D^{*+}#rightarrow#pi^{+}K^{#minus}K^{+}",                  3007},
  {DstarToPiKKPi0,      Mothers[Dstar],    "PiKKPi0",      BRsD0_PYTHIA[DstarToPiKKPi0],          BRsD0_PDG[DstarToPiKKPi0],        "D^{*+}#rightarrow#pi^{+}K^{#minus}K^{+}#pi^{0}",           3008},
  {DstarToPiPiPi,       Mothers[Dstar],    "PiPiPi",       BRsD0_PYTHIA[DstarToPiPiPi],           BRsD0_PDG[DstarToPiPiPi],         "D^{*+}#rightarrow#pi^{+}#pi^{#minus}#pi^{+}",              3009},
  {DstarToPiPiPiPi0,    Mothers[Dstar],    "PiPiPiPi0",    BRsD0_PYTHIA[DstarToPiPiPiPi0],        BRsD0_PDG[DstarToPiPiPiPi0],      "D^{*+}#rightarrow#pi^{+}#pi^{#minus}#pi^{+}#pi^{0}",       3010},

  {LcToPKPi,            Mothers[Lc],       "PKPi",         BRsLc_PYTHIA[LcToPKPi],                BRsLc_PDG[LcToPKPi],              "#Lambda_{c}^{+}#rightarrowpK^{#minus}#pi^{+}",             3001},
  {LcToPKPiPi0,         Mothers[Lc],       "PKPiPi0",      BRsLc_PYTHIA[LcToPKPiPi0],             BRsLc_PDG[LcToPKPiPi0],           "#Lambda_{c}^{+}#rightarrowpK^{#minus}#pi^{+}#pi^{0}",      3004},
  {LcToPPiPi,           Mothers[Lc],       "PPiPi",        BRsLc_PYTHIA[LcToPPiPi],               BRsLc_PDG[LcToPPiPi],             "#Lambda_{c}^{+}#rightarrowp#pi^{#minus}#pi^{+}",           3005},
  {LcToPKK,             Mothers[Lc],       "PKK",          BRsLc_PYTHIA[LcToPKK],                 BRsLc_PDG[LcToPKK],               "#Lambda_{c}^{+}#rightarrowpK^{#minus}K^{+}",               3006},

  {XicToPKPi,           Mothers[Xic],      "PKPi",         BRsXic_PYTHIA[XicToPKPi],              BRsXic_PDG[XicToPKPi],            "#Xi_{c}^{+}#rightarrowpK^{#minus}#pi^{+}",                 3004},
  {XicToPKK,            Mothers[Xic],      "PKK",          BRsXic_PYTHIA[XicToPKK],               BRsXic_PDG[XicToPKK],             "#Xi_{c}^{+}#rightarrowpK^{#minus}K^{+}",                   3005},
  {XicToSPiPi,          Mothers[Xic],      "SPiPi",        BRsXic_PYTHIA[XicToSPiPi],             BRsXic_PDG[XicToSPiPi],           "#Xi_{c}^{+}#rightarrow#Sigma^{+}#pi^{#minus}#pi^{+}",      3006},
};

//__________________________________________________________
//__________________________________________________________
//__________________________________________________________
void corrBkgLc(bool scaleByBrs = true,
               bool applyBdt = false,
               std::string filename = "/home/oleksii/alidir/sandbox/AO2D.root",
               std::vector<float> ptMins = {1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 16}, //{1, 2, 4, 6, 8},
               std::vector<float> ptMaxs = {2, 3, 4, 5, 6, 7, 8, 10, 12, 16, 24}, //{2, 4, 6, 8, 12},
               std::vector<float> bdt = {0.06, 0.04, 0.04, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.15, 0.2})
{
    const size_t nDecays{Decays.size()};
    const size_t nMothers{Mothers.size()};
    const int indexLcToPKPi = std::distance(Decays.begin(), std::find_if(Decays.begin(), Decays.end(), [](const Decay& decay) { return decay.id_ == LcToPKPi; }));
    const int indexDplusToPiKPi = std::distance(Decays.begin(), std::find_if(Decays.begin(), Decays.end(), [](const Decay& decay) { return decay.id_ == DplusToPiKPi; }));

    const std::vector<std::string> dirNames = GetDFNames(filename)/*{"DF_2263624225188351"}*/;
    TFile* file = TFile::Open(filename.c_str());

    const float minX = 2.1;
    const float maxX = 2.5;
    const int binsX = 200;

    gStyle->SetPalette(kVisibleSpectrum);
    auto cols = TColor::GetPalette();
    std::vector<int> colors = {};
    for(int i=0; i<NChannelsMain; i++) {
        colors.push_back(cols.At(i) + 9*i);
    }

    std::vector<std::vector<TH1D*>> histos;
    histos.resize(nDecays);

    std::vector<TH1D*> histos_bkgSum = {};


    int col, row, w, h;
    setRowCols(ptMins.size(), col, row, w, h);

    std::vector<TCanvas*> can;
    for(size_t iMother=0; iMother<nMothers+1; ++iMother) {
      const std::string cName = iMother == nMothers ? "can_allChannels" : "can_" + Mothers.at(iMother).name_ + "Channels";
      can.push_back(new TCanvas(cName.c_str(), cName.c_str(), w, h));
      can.back()->Divide(col, row);
    }

    /// loop over pt intervals
    for(int ptBin=0; ptBin<ptMins.size(); ptBin++) {
        const float ptMin = ptMins.at(ptBin);
        const float ptMax = ptMaxs.at(ptBin);
        std::string title = to_string_fixed_digits(ptMin, 1) + std::string(" < #it{p}_{T} < ") + to_string_fixed_digits(ptMax, 1) + std::string(" GeV/#it{c}");
//         std::vector<float> max_values = {};
        double maxY{1.};

        /// setup histograms
        /// selections
        // pt interval
        std::string cuts_base = "";
        cuts_base = cuts_base + std::to_string(ptMin) + std::string(" < fPt && fPt < ") + std::to_string(ptMax);
        // bdt
        std::string cuts_BDT = applyBdt ? std::string(" && fMlScoreFirstClass < ") + std::to_string(bdt.at(ptBin)) : std::string("");
        cuts_base = cuts_base + cuts_BDT;
        //

        std::vector<std::string> histoNames;
        std::vector<std::string> cuts;

        for(size_t iDecay=0; iDecay<nDecays; ++iDecay) {
          const std::string decayFormula = Decays.at(iDecay).mother_.name_ + "To" + Decays.at(iDecay).daughters_;
          histoNames.push_back("histos_" + decayFormula + "_ptBin" + std::to_string(ptBin));
          histos.at(iDecay).push_back(new TH1D(histoNames.at(iDecay).c_str(), histoNames.at(iDecay).c_str(), binsX, minX, maxX));

          const std::string cut = cuts_base + " && (fFlagMc == " + std::to_string(-Decays.at(iDecay).id_) + " || fFlagMc == " + std::to_string(Decays.at(iDecay).id_) + ")";
          cuts.push_back(cut);

          std::cout << cuts.back() << "\n";
        }

        std::string name_bkgSum = std::string("histos_bkgSum_ptBin") + std::to_string(ptBin);
        histos_bkgSum.push_back(new TH1D(name_bkgSum.c_str(), name_bkgSum.c_str(), binsX, minX, maxX));

        /// Open the file and get the tree
        for(const auto& dirName : dirNames) {
            std::cout << "Processing " << dirName << "\n";
            TTree* tree = file->Get<TTree>((dirName + "/O2hfcandlclite").c_str());

            std::cout << "tree->GetEntries() = " << tree->GetEntries() << "\n";

            for(size_t iDecay=0; iDecay<nDecays; ++iDecay) {
              tree->Draw(Form("fM>>+%s", histoNames.at(iDecay).c_str()), cuts.at(iDecay).c_str(), "goff");
            }
        }

        /// scale by BR values
        if (scaleByBrs) {
          for(size_t iDecay=0; iDecay<nDecays; ++iDecay) {
            const auto& decay = Decays.at(iDecay);
            histos.at(iDecay).back()->Scale(decay.br_pdg_ / (decay.mother_.pythia_br_scaling_factor_ * decay.br_pythia_));
          }
        }

        if(indexLcToPKPi < nDecays) std::cout << "Lc: " << histos.at(indexLcToPKPi).back()->Integral() << std::endl;
        if(indexDplusToPiKPi < nDecays) std::cout << "D+: " << histos.at(indexDplusToPiKPi).back()->Integral() << std::endl;

        /// create the template with the sum of correlated background sources
        for(size_t iDecay=0; iDecay<nDecays; ++iDecay) {
          const auto& decay = Decays.at(iDecay);
          const auto& histo = histos.at(iDecay).back();
          if(decay.id_ != LcToPKPi) histos_bkgSum.back()->Add(histo); // this is the signal!

          const auto& color = colors.at(decay.id_-1);
          histo->SetTitle(decay.formula_.c_str());
          histo->SetLineColor(color);
          histo->SetLineWidth(2);
          histo->SetMarkerColor(color);
          histo->SetFillColor(color);
          histo->SetFillStyle(decay.fill_style_);
          maxY = std::max(maxY, histo->GetMaximum());
        }
        const float binwidth = histos.at(0).back()->GetBinWidth(1) * 1000; // MeV/c2


        // sum
        histos_bkgSum.back()->SetTitle("sum of backgrounds");
        histos_bkgSum.back()->SetLineColor(kGray+1);
        histos_bkgSum.back()->SetLineWidth(2);
        histos_bkgSum.back()->SetMarkerColor(histos_bkgSum.back()->GetLineColor());

        /// Draw stuff
        // D+ channels
        std::vector<TH1F*> frames;
        frames.resize(nMothers+1);
        for(int iMother=0; iMother<nMothers+1; ++iMother) {
          auto& frame = frames.at(iMother);
          frame = can.at(iMother)->cd(ptBin+1)->DrawFrame(2.1, 0.00000001, 2.5, maxY, title.c_str());
          frame->GetXaxis()->SetTitle("#it{M}(pK#pi) (GeV/#it{c}^{2})");
          frame->GetYaxis()->SetTitle(Form("counts/%.0f MeV/#it{c}^{2}", binwidth));
          frame->SetLineWidth(0);
          for(int iDecay=0; iDecay<nDecays; ++iDecay) {
            const auto& decay = Decays.at(iDecay);
            if(iMother != nMothers && decay.mother_.id_ != Mothers.at(iMother).id_) continue;
            histos.at(iDecay).back()->Draw("samehist");
          }
          if(iMother == nMothers) histos_bkgSum.back()->Draw("samehist");
          gPad->Update();
          gPad->SetBottomMargin(0.15);
          gPad->SetLeftMargin(0.15);
          gPad->SetTicks();
          gPad->Update();
        }

        if(ptBin == row*col-2 || ptBin == row*col-1) {
            int ptBin_legend = ptBin;
            float minY_legend = 0.6;
            float maxY_legend = 0.8;
            if (ptBin == row*col-2) {
                ptBin_legend++;
                minY_legend = 0.2;
                maxY_legend = 0.85;
            }
            std::vector<TLegend*> legends;
            legends.resize(nMothers+1);
            for(int iMother=0; iMother<nMothers+1; ++iMother) {
              auto& leg = legends.at(iMother);
              leg = new TLegend(0.2, minY_legend, 0.8, maxY_legend);
              leg->SetFillStyle(0);
              leg->SetBorderSize(0);
              for(int iDecay=0; iDecay<nDecays; ++iDecay) {
                const auto& decay = Decays.at(iDecay);
                if(iMother != nMothers && decay.mother_.id_ != Mothers.at(iMother).id_) continue;
                leg->AddEntry(histos.at(iDecay).back());
              }
              if(iMother == nMothers) {
                leg->SetNColumns(2);
                leg->AddEntry(histos_bkgSum.back());
              }
              can.at(iMother)->cd(ptBin_legend+1);
              leg->Draw();
              gPad->Update();
            }
        }
        

    } /// end loop over pt intervals

    /// save stuff
    TFile* fout = new TFile(Form("file_corrBkgPKPi_MC%s%s.root", scaleByBrs ? "_brScaled" : "", applyBdt ? "_bdtCuts" : ""), "recreate");
    for(int iMother=0; iMother<nMothers+1; ++iMother) {
      can.at(iMother)->Write();
    }
    for(int ptBin=0; ptBin<ptMins.size(); ptBin++) {
        fout ->cd();
        TDirectory* td = fout->mkdir(Form("ptBin%d", ptBin));
        td->cd();
        for(int iDecay=0; iDecay<nDecays; ++iDecay) {
          histos.at(iDecay).at(ptBin)->Write();
        }
        histos_bkgSum.at(ptBin)->Write();
    }
    fout->Close();

    file->Close();

    return;
}
//__________________________________________________________
//__________________________________________________________
//__________________________________________________________
void setRowCols(const int size, int& rows, int& cols, int& w, int& h) {
    if(size == 1) {
        rows = 1;
        cols = 1;
        w = 800;
        h = 800;
    } else if(size == 2) {
        rows = 2;
        cols = 1;
        w = 1400;
        h = 700;
    } else if(size <= 4) {
        rows = 2;
        cols = 2;
        w = 1000;
        h = 1000;
    } else if(size <= 6) {
        rows = 3;
        cols = 2;
        w = 1500;
        h = 1000;
    } else if(size <= 8) {
        rows = 2;
        cols = 4;
         w = 1600;
         h = 800;
    } else if(size <= 12) {
        rows = 4;
        cols = 3;
        w = 1600;
        h = 1200;
    } else if(size <= 20) {
        rows = 4;
        cols = 5;
        w = 2000;
        h = 1600;
    } else if(size <= 25) {
        rows = 5;
        cols = 5;
        w = 1200;
        h = 800;
    } else {
        rows = 6;
        cols = 6;
        w = 1200;
        h = 1200;
    }
    
    return;
}
//__________________________________________________________
//__________________________________________________________
//__________________________________________________________
template <typename I>
std::string to_string_fixed_digits(I value,int digitsCount)
{
    std::stringstream stream;
    stream << std::fixed << std::setprecision(digitsCount) <<value;
    return stream.str();
}
//__________________________________________________________
//__________________________________________________________
//__________________________________________________________
std::vector<std::string> GetDFNames(const std::string& fileName) {
  TFile* fileIn = TFile::Open(fileName.c_str());
  if(fileIn == nullptr) {
    throw std::runtime_error("fileIn == nullptr");
  }

  std::vector<std::string> result;
  auto lok = fileIn->GetListOfKeys();
  for(const auto& k : *lok) {
    const std::string dirname = k->GetName();
    if(dirname.substr(0, 2) != "DF") continue;
    result.emplace_back(dirname);
  }
  fileIn->Close();

  return result;
}
