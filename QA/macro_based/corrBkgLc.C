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

    const std::vector<std::string> dirNames = GetDFNames(filename);
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

    std::vector<TH1D*> histos_DplusToPiKPi = {};
    std::vector<TH1D*> histos_DplusToPiKPiPi0 = {};
    std::vector<TH1D*> histos_DplusToPiPiPi = {};
    std::vector<TH1D*> histos_DplusToPiKK = {};
    std::vector<TH1D*> histos_DsToPiKK = {};
    std::vector<TH1D*> histos_DsToPiKKPi0 = {};
    std::vector<TH1D*> histos_DsToPiPiK = {};
    std::vector<TH1D*> histos_DsToPiPiPi = {};
    std::vector<TH1D*> histos_DsToPiPiPiPi0 = {};
    std::vector<TH1D*> histos_DstarToPiKPi = {};
    std::vector<TH1D*> histos_DstarToPiKPiPi0 = {};
    std::vector<TH1D*> histos_DstarToPiKPiPi0Pi0 = {};
    std::vector<TH1D*> histos_DstarToPiKK = {};
    std::vector<TH1D*> histos_DstarToPiKKPi0 = {};
    std::vector<TH1D*> histos_DstarToPiPiPi = {};
    std::vector<TH1D*> histos_DstarToPiPiPiPi0 = {};
    std::vector<TH1D*> histos_LcToPKPi = {};
    std::vector<TH1D*> histos_LcToPKPiPi0 = {};
    std::vector<TH1D*> histos_LcToPPiPi = {};
    std::vector<TH1D*> histos_LcToPKK = {};
    std::vector<TH1D*> histos_XicToPKPi = {};
    std::vector<TH1D*> histos_XicToPKK = {};
    std::vector<TH1D*> histos_XicToSPiPi = {};
    std::vector<TH1D*> histos_bkgSum = {};


    int col, row, w, h;
    setRowCols(ptMins.size(), col, row, w, h);
    TCanvas* can_DplusChannels = new TCanvas("can_DplusChannels", "can_DplusChannels", w, h);
    can_DplusChannels->Divide(col, row);
    TCanvas* can_DsChannels = new TCanvas("can_DsChannels", "can_DsChannels", w, h);
    can_DsChannels->Divide(col, row);
    TCanvas* can_DstarChannels = new TCanvas("can_DstarChannels", "can_DstarChannels", w, h);
    can_DstarChannels->Divide(col, row);
    TCanvas* can_LcChannels = new TCanvas("can_LcChannels", "can_LcChannels", w, h);
    can_LcChannels->Divide(col, row);
    TCanvas* can_XicChannels = new TCanvas("can_XicChannels", "can_XicChannels", w, h);
    can_XicChannels->Divide(col, row);
    TCanvas* can_allChannels = new TCanvas("can_allChannels", "can_allChannels", w, h);
    can_allChannels->Divide(col, row);
    float maxY = 3000;

    /// loop over pt intervals
    for(int ptBin=0; ptBin<ptMins.size(); ptBin++) {
        const float ptMin = ptMins.at(ptBin);
        const float ptMax = ptMaxs.at(ptBin);
        std::string title = to_string_fixed_digits(ptMin, 1) + std::string(" < #it{p}_{T} < ") + to_string_fixed_digits(ptMax, 1) + std::string(" GeV/#it{c}");
        std::vector<float> max_values = {};

        /// setup histograms
        std::string name_DplusToPiKPi = std::string("histos_DplusToPiKPi_ptBin") + std::to_string(ptBin);
        std::string name_DplusToPiKPiPi0 = std::string("histos_DplusToPiKPiPi0_ptBin") + std::to_string(ptBin);
        std::string name_DplusToPiPiPi = std::string("histos_DplusToPiPiPi_ptBin") + std::to_string(ptBin);
        std::string name_DplusToPiKK = std::string("histos_DplusToPiKK_ptBin") + std::to_string(ptBin);
        std::string name_DsToPiKK = std::string("histos_DsToPiKK_ptBin") + std::to_string(ptBin);
        std::string name_DsToPiKKPi0 = std::string("histos_DsToPiKKPi0_ptBin") + std::to_string(ptBin);
        std::string name_DsToPiPiK = std::string("histos_DsToPiPiK_ptBin") + std::to_string(ptBin);
        std::string name_DsToPiPiPi = std::string("histos_DsToPiPiPi_ptBin") + std::to_string(ptBin);
        std::string name_DsToPiPiPiPi0 = std::string("histos_DsToPiPiPiPi0_ptBin") + std::to_string(ptBin);
        std::string name_DstarToPiKPi = std::string("histos_DstarToPiKPi_ptBin") + std::to_string(ptBin);
        std::string name_DstarToPiKPiPi0 = std::string("histos_DstarToPiKPiPi0_ptBin") + std::to_string(ptBin);
        std::string name_DstarToPiKPiPi0Pi0 = std::string("histos_DstarToPiKPiPi0Pi0_ptBin") + std::to_string(ptBin);
        std::string name_DstarToPiKK = std::string("histos_DstarToPiKK_ptBin") + std::to_string(ptBin);
        std::string name_DstarToPiKKPi0 = std::string("histos_DstarToPiKKPi0_ptBin") + std::to_string(ptBin);
        std::string name_DstarToPiPiPi = std::string("histos_DstarToPiPiPi_ptBin") + std::to_string(ptBin);
        std::string name_DstarToPiPiPiPi0 = std::string("histos_DstarToPiPiPiPi0_ptBin") + std::to_string(ptBin);
        std::string name_LcToPKPi = std::string("histos_LcToPKPi_ptBin") + std::to_string(ptBin);
        std::string name_LcToPKPiPi0 = std::string("histos_LcToPKPiPi0_ptBin") + std::to_string(ptBin);
        std::string name_LcToPPiPi = std::string("histos_LcToPPiPi_ptBin") + std::to_string(ptBin);
        std::string name_LcToPKK = std::string("histos_LcToPKK_ptBin") + std::to_string(ptBin);
        std::string name_XicToPKPi = std::string("histos_XicToPKPi_ptBin") + std::to_string(ptBin);
        std::string name_XicToPKK = std::string("histos_XicToPKK_ptBin") + std::to_string(ptBin);
        std::string name_XicToSPiPi = std::string("histos_XicToSPiPi_ptBin") + std::to_string(ptBin);
        std::string name_bkgSum = std::string("histos_bkgSum_ptBin") + std::to_string(ptBin);
        histos_DplusToPiKPi.push_back(new TH1D(name_DplusToPiKPi.c_str(), name_DplusToPiKPi.c_str(), binsX, minX, maxX));
        histos_DplusToPiKPiPi0.push_back(new TH1D(name_DplusToPiKPiPi0.c_str(), name_DplusToPiKPiPi0.c_str(), binsX, minX, maxX));
        histos_DplusToPiPiPi.push_back(new TH1D(name_DplusToPiPiPi.c_str(), name_DplusToPiPiPi.c_str(), binsX, minX, maxX));
        histos_DplusToPiKK.push_back(new TH1D(name_DplusToPiKK.c_str(), name_DplusToPiKK.c_str(), binsX, minX, maxX));
        histos_DsToPiKK.push_back(new TH1D(name_DsToPiKK.c_str(), name_DsToPiKK.c_str(), binsX, minX, maxX));
        histos_DsToPiKKPi0.push_back(new TH1D(name_DsToPiKKPi0.c_str(), name_DsToPiKKPi0.c_str(), binsX, minX, maxX));
        histos_DsToPiPiK.push_back(new TH1D(name_DsToPiPiK.c_str(), name_DsToPiPiK.c_str(), binsX, minX, maxX));
        histos_DsToPiPiPi.push_back(new TH1D(name_DsToPiPiPi.c_str(), name_DsToPiPiPi.c_str(), binsX, minX, maxX));
        histos_DsToPiPiPiPi0.push_back(new TH1D(name_DsToPiPiPiPi0.c_str(), name_DsToPiPiPiPi0.c_str(), binsX, minX, maxX));
        histos_DstarToPiKPi.push_back(new TH1D(name_DstarToPiKPi.c_str(), name_DstarToPiKPi.c_str(), binsX, minX, maxX));
        histos_DstarToPiKPiPi0.push_back(new TH1D(name_DstarToPiKPiPi0.c_str(), name_DstarToPiKPiPi0.c_str(), binsX, minX, maxX));
        histos_DstarToPiKPiPi0Pi0.push_back(new TH1D(name_DstarToPiKPiPi0Pi0.c_str(), name_DstarToPiKPiPi0Pi0.c_str(), binsX, minX, maxX));
        histos_DstarToPiKK.push_back(new TH1D(name_DstarToPiKK.c_str(), name_DstarToPiKK.c_str(), binsX, minX, maxX));
        histos_DstarToPiKKPi0.push_back(new TH1D(name_DstarToPiKKPi0.c_str(), name_DstarToPiKKPi0.c_str(), binsX, minX, maxX));
        histos_DstarToPiPiPi.push_back(new TH1D(name_DstarToPiPiPi.c_str(), name_DstarToPiPiPi.c_str(), binsX, minX, maxX));
        histos_DstarToPiPiPiPi0.push_back(new TH1D(name_DstarToPiPiPiPi0.c_str(), name_DstarToPiPiPiPi0.c_str(), binsX, minX, maxX));
        histos_LcToPKPi.push_back(new TH1D(name_LcToPKPi.c_str(), name_LcToPKPi.c_str(), binsX, minX, maxX));
        histos_LcToPKPiPi0.push_back(new TH1D(name_LcToPKPiPi0.c_str(), name_LcToPKPiPi0.c_str(), binsX, minX, maxX));
        histos_LcToPPiPi.push_back(new TH1D(name_LcToPPiPi.c_str(), name_LcToPPiPi.c_str(), binsX, minX, maxX));
        histos_LcToPKK.push_back(new TH1D(name_LcToPKK.c_str(), name_LcToPKK.c_str(), binsX, minX, maxX));
        histos_XicToPKPi.push_back(new TH1D(name_XicToPKPi.c_str(), name_XicToPKPi.c_str(), binsX, minX, maxX));
        histos_XicToPKK.push_back(new TH1D(name_XicToPKK.c_str(), name_XicToPKK.c_str(), binsX, minX, maxX));
        histos_XicToSPiPi.push_back(new TH1D(name_XicToSPiPi.c_str(), name_XicToSPiPi.c_str(), binsX, minX, maxX));
        histos_bkgSum.push_back(new TH1D(name_bkgSum.c_str(), name_bkgSum.c_str(), binsX, minX, maxX));

        /// selections
        // pt interval
        std::string cuts_base = "";
        cuts_base = cuts_base + std::to_string(ptMin) + std::string(" < fPt && fPt < ") + std::to_string(ptMax);
        // bdt
        std::string cuts_BDT = applyBdt ? std::string(" && fMlScoreFirstClass < ") + std::to_string(bdt.at(ptBin)) : std::string("");
        cuts_base = cuts_base + cuts_BDT;
        //
        std::string cuts_DplusToPiKPi = cuts_base + std::string(" && (fFlagMc == ") + std::to_string(-DplusToPiKPi) + std::string(" || fFlagMc == ") + std::to_string(DplusToPiKPi) + std::string(")");
        std::string cuts_DplusToPiKPiPi0 = cuts_base + std::string(" && (fFlagMc == ") + std::to_string(-DplusToPiKPiPi0) + std::string(" || fFlagMc == ") + std::to_string(DplusToPiKPiPi0) + std::string(")");
        std::string cuts_DplusToPiPiPi = cuts_base + std::string(" && (fFlagMc == ") + std::to_string(-DplusToPiPiPi) + std::string(" || fFlagMc == ") + std::to_string(DplusToPiPiPi) + std::string(")");
        std::string cuts_DplusToPiKK = cuts_base + std::string(" && (fFlagMc == ") + std::to_string(-DplusToPiKK) + std::string(" || fFlagMc == ") + std::to_string(DplusToPiKK) + std::string(")");
        std::string cuts_DsToPiKK = cuts_base + std::string(" && (fFlagMc == ") + std::to_string(-DsToPiKK) + std::string(" || fFlagMc == ") + std::to_string(DsToPiKK) + std::string(")");
        std::string cuts_DsToPiKKPi0 = cuts_base + std::string(" && (fFlagMc == ") + std::to_string(-DsToPiKKPi0) + std::string(" || fFlagMc == ") + std::to_string(DsToPiKKPi0) + std::string(")");
        std::string cuts_DsToPiPiK = cuts_base + std::string(" && (fFlagMc == ") + std::to_string(-DsToPiPiK) + std::string(" || fFlagMc == ") + std::to_string(DsToPiPiK) + std::string(")");
        std::string cuts_DsToPiPiPi = cuts_base + std::string(" && (fFlagMc == ") + std::to_string(-DsToPiPiPi) + std::string(" || fFlagMc == ") + std::to_string(DsToPiPiPi) + std::string(")");
        std::string cuts_DsToPiPiPiPi0 = cuts_base + std::string(" && (fFlagMc == ") + std::to_string(-DsToPiPiPiPi0) + std::string(" || fFlagMc == ") + std::to_string(DsToPiPiPiPi0) + std::string(")");
        std::string cuts_DstarToPiKPi = cuts_base + std::string(" && (fFlagMc == ") + std::to_string(-DstarToPiKPi) + std::string(" || fFlagMc == ") + std::to_string(DstarToPiKPi) + std::string(")");
        std::string cuts_DstarToPiKPiPi0 = cuts_base + std::string(" && (fFlagMc == ") + std::to_string(-DstarToPiKPiPi0) + std::string(" || fFlagMc == ") + std::to_string(DstarToPiKPiPi0) + std::string(")");
        std::string cuts_DstarToPiKPiPi0Pi0 = cuts_base + std::string(" && (fFlagMc == ") + std::to_string(-DstarToPiKPiPi0Pi0) + std::string(" || fFlagMc == ") + std::to_string(DstarToPiKPiPi0Pi0) + std::string(")");
        std::string cuts_DstarToPiKK = cuts_base + std::string(" && (fFlagMc == ") + std::to_string(-DstarToPiKK) + std::string(" || fFlagMc == ") + std::to_string(DstarToPiKK) + std::string(")");
        std::string cuts_DstarToPiKKPi0 = cuts_base + std::string(" && (fFlagMc == ") + std::to_string(-DstarToPiKKPi0) + std::string(" || fFlagMc == ") + std::to_string(DstarToPiKKPi0) + std::string(")");
        std::string cuts_DstarToPiPiPi = cuts_base + std::string(" && (fFlagMc == ") + std::to_string(-DstarToPiPiPi) + std::string(" || fFlagMc == ") + std::to_string(DstarToPiPiPi) + std::string(")");
        std::string cuts_DstarToPiPiPiPi0 = cuts_base + std::string(" && (fFlagMc == ") + std::to_string(-DstarToPiPiPiPi0) + std::string(" || fFlagMc == ") + std::to_string(DstarToPiPiPiPi0) + std::string(")");
        std::string cuts_LcToPKPi = cuts_base + std::string(" && (fFlagMc == ") + std::to_string(-LcToPKPi) + std::string(" || fFlagMc == ") + std::to_string(LcToPKPi) + std::string(")");
        std::string cuts_LcToPKPiPi0 = cuts_base + std::string(" && (fFlagMc == ") + std::to_string(-LcToPKPiPi0) + std::string(" || fFlagMc == ") + std::to_string(LcToPKPiPi0) + std::string(")");
        std::string cuts_LcToPPiPi = cuts_base + std::string(" && (fFlagMc == ") + std::to_string(-LcToPPiPi) + std::string(" || fFlagMc == ") + std::to_string(LcToPPiPi) + std::string(")");
        std::string cuts_LcToPKK = cuts_base + std::string(" && (fFlagMc == ") + std::to_string(-LcToPKK) + std::string(" || fFlagMc == ") + std::to_string(LcToPKK) + std::string(")");
        std::string cuts_XicToPKPi = cuts_base + std::string(" && (fFlagMc == ") + std::to_string(-XicToPKPi) + std::string(" || fFlagMc == ") + std::to_string(XicToPKPi) + std::string(")");
        std::string cuts_XicToPKK = cuts_base + std::string(" && (fFlagMc == ") + std::to_string(-XicToPKK) + std::string(" || fFlagMc == ") + std::to_string(XicToPKK) + std::string(")");
        std::string cuts_XicToSPiPi = cuts_base + std::string(" && (fFlagMc == ") + std::to_string(-XicToSPiPi) + std::string(" || fFlagMc == ") + std::to_string(XicToSPiPi) + std::string(")");

        std::cout << std::endl;
        std::cout << cuts_DplusToPiKPi << std::endl;
        std::cout << cuts_DplusToPiKPiPi0 << std::endl;
        std::cout << cuts_DplusToPiPiPi << std::endl;
        std::cout << cuts_DplusToPiKK << std::endl;
        std::cout << cuts_DsToPiKK << std::endl;
        std::cout << cuts_DsToPiKKPi0 << std::endl;
        std::cout << cuts_DsToPiPiK << std::endl;
        std::cout << cuts_DsToPiPiPi << std::endl;
        std::cout << cuts_DsToPiPiPiPi0 << std::endl;
        std::cout << cuts_DstarToPiKPi << std::endl;
        std::cout << cuts_DstarToPiKPiPi0 << std::endl;
        std::cout << cuts_DstarToPiKPiPi0Pi0 << std::endl;
        std::cout << cuts_DstarToPiKK << std::endl;
        std::cout << cuts_DstarToPiKKPi0 << std::endl;
        std::cout << cuts_DstarToPiPiPi << std::endl;
        std::cout << cuts_DstarToPiPiPiPi0 << std::endl;
        std::cout << cuts_LcToPKPi << std::endl;
        std::cout << cuts_LcToPKPiPi0 << std::endl;
        std::cout << cuts_LcToPPiPi << std::endl;
        std::cout << cuts_LcToPKK << std::endl;
        std::cout << cuts_XicToPKPi << std::endl;
        std::cout << cuts_XicToPKK << std::endl;
        std::cout << cuts_XicToSPiPi << std::endl;



        /// Open the file and get the tree
        for(const auto& dirName : dirNames) {
            std::cout << "Processing " << dirName << "\n";
            TTree* tree = file->Get<TTree>((dirName + "/O2hfcandlclite").c_str());

            std::cout << "tree->GetEntries() = " << tree->GetEntries() << "\n";

            /// TTree projections
            // D+ channels
            tree->Draw(Form("fM>>%s", name_DplusToPiKPi.c_str()), cuts_DplusToPiKPi.c_str(), "goff");
            tree->Draw(Form("fM>>%s", name_DplusToPiKPiPi0.c_str()), cuts_DplusToPiKPiPi0.c_str(), "goff");
            tree->Draw(Form("fM>>%s", name_DplusToPiPiPi.c_str()), cuts_DplusToPiPiPi.c_str(), "goff");
            tree->Draw(Form("fM>>%s", name_DplusToPiKK.c_str()), cuts_DplusToPiKK.c_str(), "goff");
            // Ds channels
            tree->Draw(Form("fM>>%s", name_DsToPiKK.c_str()), cuts_DsToPiKK.c_str(), "goff");
            tree->Draw(Form("fM>>%s", name_DsToPiKKPi0.c_str()), cuts_DsToPiKKPi0.c_str(), "goff");
            tree->Draw(Form("fM>>%s", name_DsToPiPiK.c_str()), cuts_DsToPiPiK.c_str(), "goff");
            tree->Draw(Form("fM>>%s", name_DsToPiPiPi.c_str()), cuts_DsToPiPiPi.c_str(), "goff");
            tree->Draw(Form("fM>>%s", name_DsToPiPiPiPi0.c_str()), cuts_DsToPiPiPiPi0.c_str(), "goff");
            // D* channels
            tree->Draw(Form("fM>>%s", name_DstarToPiKPi.c_str()), cuts_DstarToPiKPi.c_str(), "goff");
            tree->Draw(Form("fM>>%s", name_DstarToPiKPiPi0.c_str()), cuts_DstarToPiKPiPi0.c_str(), "goff");
            tree->Draw(Form("fM>>%s", name_DstarToPiKPiPi0Pi0.c_str()), cuts_DstarToPiKPiPi0Pi0.c_str(), "goff");
            tree->Draw(Form("fM>>%s", name_DstarToPiKK.c_str()), cuts_DstarToPiKK.c_str(), "goff");
            tree->Draw(Form("fM>>%s", name_DstarToPiKKPi0.c_str()), cuts_DstarToPiKKPi0.c_str(), "goff");
            tree->Draw(Form("fM>>%s", name_DstarToPiPiPi.c_str()), cuts_DstarToPiPiPi.c_str(), "goff");
            tree->Draw(Form("fM>>%s", name_DstarToPiPiPiPi0.c_str()), cuts_DstarToPiPiPiPi0.c_str(), "goff");
            // Lc channels
            tree->Draw(Form("fM>>%s", name_LcToPKPi.c_str()), cuts_LcToPKPi.c_str(), "goff");
            tree->Draw(Form("fM>>%s", name_LcToPKPiPi0.c_str()), cuts_LcToPKPiPi0.c_str(), "goff");
            tree->Draw(Form("fM>>%s", name_LcToPPiPi.c_str()), cuts_LcToPPiPi.c_str(), "goff");
            tree->Draw(Form("fM>>%s", name_LcToPKK.c_str()), cuts_LcToPKK.c_str(), "goff");
            // Xic channels
            tree->Draw(Form("fM>>%s", name_XicToPKPi.c_str()), cuts_XicToPKPi.c_str(), "goff");
            tree->Draw(Form("fM>>%s", name_XicToPKK.c_str()), cuts_XicToPKK.c_str(), "goff");
            tree->Draw(Form("fM>>%s", name_XicToSPiPi.c_str()), cuts_XicToSPiPi.c_str(), "goff");

        }


        /// scale by BR values
        if (scaleByBrs) {
            // D+ channels
            histos_DplusToPiKPi.back()->Scale(      BRsDplus_PDG[DplusToPiKPi]      / (pythiaBRScalingFactorDPlus * BRsDplus_PYTHIA[DplusToPiKPi]));
            histos_DplusToPiKPiPi0.back()->Scale(   BRsDplus_PDG[DplusToPiKPiPi0]   / (pythiaBRScalingFactorDPlus * BRsDplus_PYTHIA[DplusToPiKPiPi0]));
            histos_DplusToPiPiPi.back()->Scale(     BRsDplus_PDG[DplusToPiPiPi]     / (pythiaBRScalingFactorDPlus * BRsDplus_PYTHIA[DplusToPiPiPi]));
            histos_DplusToPiKK.back()->Scale(       BRsDplus_PDG[DplusToPiKK]       / (pythiaBRScalingFactorDPlus * BRsDplus_PYTHIA[DplusToPiKK]));
            // Ds channels
            histos_DsToPiKK.back()->Scale(      BRsDs_PDG[DsToPiKK]      / (pythiaBRScalingFactorDs * BRsDs_PYTHIA[DsToPiKK]));
            histos_DsToPiKKPi0.back()->Scale(   BRsDs_PDG[DsToPiKKPi0]   / (pythiaBRScalingFactorDs * BRsDs_PYTHIA[DsToPiKKPi0]));
            histos_DsToPiPiK.back()->Scale(     BRsDs_PDG[DsToPiPiK]     / (pythiaBRScalingFactorDs * BRsDs_PYTHIA[DsToPiPiK]));
            histos_DsToPiPiPi.back()->Scale(    BRsDs_PDG[DsToPiPiPi]    / (pythiaBRScalingFactorDs * BRsDs_PYTHIA[DsToPiPiPi]));
            histos_DsToPiPiPiPi0.back()->Scale( BRsDs_PDG[DsToPiPiPiPi0] / (pythiaBRScalingFactorDs * BRsDs_PYTHIA[DsToPiPiPiPi0]));
            // Dstar channels
            histos_DstarToPiKPi.back()->Scale(          BRsD0_PDG[DstarToPiKPi]         / (pythiaBRScalingFactorDstar * BRsD0_PYTHIA[DstarToPiKPi]));
            histos_DstarToPiKPiPi0.back()->Scale(       BRsD0_PDG[DstarToPiKPiPi0]      / (pythiaBRScalingFactorDstar * BRsD0_PYTHIA[DstarToPiKPiPi0]));
            histos_DstarToPiKPiPi0Pi0.back()->Scale(    BRsD0_PDG[DstarToPiKPiPi0Pi0]   / (pythiaBRScalingFactorDstar * BRsD0_PYTHIA[DstarToPiKPiPi0Pi0]));
            histos_DstarToPiKK.back()->Scale(           BRsD0_PDG[DstarToPiKK]          / (pythiaBRScalingFactorDstar * BRsD0_PYTHIA[DstarToPiKK]));
            histos_DstarToPiKKPi0.back()->Scale(        BRsD0_PDG[DstarToPiKKPi0]       / (pythiaBRScalingFactorDstar * BRsD0_PYTHIA[DstarToPiKKPi0]));
            histos_DstarToPiPiPi.back()->Scale(         BRsD0_PDG[DstarToPiPiPi]        / (pythiaBRScalingFactorDstar * BRsD0_PYTHIA[DstarToPiPiPi]));
            histos_DstarToPiPiPiPi0.back()->Scale(      BRsD0_PDG[DstarToPiPiPiPi0]     / (pythiaBRScalingFactorDstar * BRsD0_PYTHIA[DstarToPiPiPiPi0]));
            // Lc channels
            histos_LcToPKPi.back()->Scale(      BRsLc_PDG[LcToPKPi]     / (pythiaBRScalingFactorLc * BRsLc_PYTHIA[LcToPKPi]));
            histos_LcToPKPiPi0.back()->Scale(   BRsLc_PDG[LcToPKPiPi0]  / (pythiaBRScalingFactorLc * BRsLc_PYTHIA[LcToPKPiPi0]));
            histos_LcToPPiPi.back()->Scale(     BRsLc_PDG[LcToPPiPi]    / (pythiaBRScalingFactorLc * BRsLc_PYTHIA[LcToPPiPi]));
            histos_LcToPKK.back()->Scale(       BRsLc_PDG[LcToPKK]      / (pythiaBRScalingFactorLc * BRsLc_PYTHIA[LcToPKK]));
            // Xic channels
            histos_XicToPKPi.back()->Scale(    BRsXic_PDG[XicToPKPi]    / (pythiaBRScalingFactorXic * BRsXic_PYTHIA[XicToPKPi]));
            histos_XicToPKK.back()->Scale(     BRsXic_PDG[XicToPKK]     / (pythiaBRScalingFactorXic * BRsXic_PYTHIA[XicToPKK]));
            histos_XicToSPiPi.back()->Scale(   BRsXic_PDG[XicToSPiPi]   / (pythiaBRScalingFactorXic * BRsXic_PYTHIA[XicToSPiPi]));
        }

        std::cout << "Lc: " << histos_LcToPKPi.back()->Integral() << std::endl;
        std::cout << "D+: " << histos_DplusToPiKPi.back()->Integral() << std::endl;

        /// create the template with the sum of correlated background sources
        histos_bkgSum.back()->Add(histos_DplusToPiKPi.back());
        histos_bkgSum.back()->Add(histos_DplusToPiKPiPi0.back());
        histos_bkgSum.back()->Add(histos_DplusToPiPiPi.back());
        histos_bkgSum.back()->Add(histos_DplusToPiKK.back());
        histos_bkgSum.back()->Add(histos_DsToPiKK.back());
        histos_bkgSum.back()->Add(histos_DsToPiKKPi0.back());
        histos_bkgSum.back()->Add(histos_DsToPiPiK.back());
        histos_bkgSum.back()->Add(histos_DsToPiPiPi.back());
        histos_bkgSum.back()->Add(histos_DsToPiPiPiPi0.back());
        histos_bkgSum.back()->Add(histos_DstarToPiKPi.back());
        histos_bkgSum.back()->Add(histos_DstarToPiKPiPi0.back());
        histos_bkgSum.back()->Add(histos_DstarToPiKPiPi0Pi0.back());
        histos_bkgSum.back()->Add(histos_DstarToPiKK.back());
        histos_bkgSum.back()->Add(histos_DstarToPiKKPi0.back());
        histos_bkgSum.back()->Add(histos_DstarToPiPiPi.back());
        histos_bkgSum.back()->Add(histos_DstarToPiPiPiPi0.back());
        //histos_bkgSum.back()->Add(histos_LcToPKPi.back());    // this is the signal!
        histos_bkgSum.back()->Add(histos_LcToPKPiPi0.back());
        histos_bkgSum.back()->Add(histos_LcToPPiPi.back());
        histos_bkgSum.back()->Add(histos_LcToPKK.back());
        histos_bkgSum.back()->Add(histos_XicToPKPi.back());
        histos_bkgSum.back()->Add(histos_XicToPKK.back());
        histos_bkgSum.back()->Add(histos_XicToSPiPi.back());

        /// set style
        // D+ channels
        histos_DplusToPiKPi.back()->SetTitle("D^{+}#rightarrowK^{#minus}#pi^{+}#pi^{+}");
        histos_DplusToPiKPi.back()->SetLineColor(colors.at(DplusToPiKPi-1));
        histos_DplusToPiKPi.back()->SetLineWidth(2);
        histos_DplusToPiKPi.back()->SetMarkerColor(histos_DplusToPiKPi.back()->GetLineColor());
        histos_DplusToPiKPi.back()->SetFillColor(histos_DplusToPiKPi.back()->GetLineColor());
        histos_DplusToPiKPi.back()->SetFillStyle(3004);
        max_values.push_back(histos_DplusToPiKPi.back()->GetMaximum());
        float binwidth = histos_DplusToPiKPi.back()->GetBinWidth(1) * 1000; // MeV/c2
        //
        histos_DplusToPiKPiPi0.back()->SetTitle("D^{+}#rightarrowK^{#minus}#pi^{+}#pi^{+}#pi^{0}");
        histos_DplusToPiKPiPi0.back()->SetLineColor(colors.at(DplusToPiKPiPi0-1));
        histos_DplusToPiKPiPi0.back()->SetLineWidth(2);
        histos_DplusToPiKPiPi0.back()->SetMarkerColor(histos_DplusToPiKPiPi0.back()->GetLineColor());
        histos_DplusToPiKPiPi0.back()->SetFillColor(histos_DplusToPiKPiPi0.back()->GetLineColor());
        histos_DplusToPiKPiPi0.back()->SetFillStyle(3005);
        max_values.push_back(histos_DplusToPiKPiPi0.back()->GetMaximum());
        //
        histos_DplusToPiPiPi.back()->SetTitle("D^{+}#rightarrow#pi^{#minus}#pi^{+}#pi^{+}");
        histos_DplusToPiPiPi.back()->SetLineColor(colors.at(DplusToPiPiPi-1));
        histos_DplusToPiPiPi.back()->SetLineWidth(2);
        histos_DplusToPiPiPi.back()->SetMarkerColor(histos_DplusToPiPiPi.back()->GetLineColor());
        histos_DplusToPiPiPi.back()->SetFillColor(histos_DplusToPiPiPi.back()->GetLineColor());
        histos_DplusToPiPiPi.back()->SetFillStyle(3006);
        max_values.push_back(histos_DplusToPiPiPi.back()->GetMaximum());
        //
        histos_DplusToPiKK.back()->SetTitle("D^{+}#rightarrow#pi^{+}K^{#minus}K^{+}");
        histos_DplusToPiKK.back()->SetLineColor(colors.at(DplusToPiKK-1));
        histos_DplusToPiKK.back()->SetLineWidth(2);
        histos_DplusToPiKK.back()->SetMarkerColor(histos_DplusToPiKK.back()->GetLineColor());
        histos_DplusToPiKK.back()->SetFillColor(histos_DplusToPiKK.back()->GetLineColor());
        histos_DplusToPiKK.back()->SetFillStyle(3007);
        max_values.push_back(histos_DplusToPiKK.back()->GetMaximum());
        // Ds channels
        histos_DsToPiKK.back()->SetTitle("D_{s}^{+}#rightarrowK^{#minus}#pi^{+}#pi^{+}");
        histos_DsToPiKK.back()->SetLineColor(colors.at(DsToPiKK-1));
        histos_DsToPiKK.back()->SetLineWidth(2);
        histos_DsToPiKK.back()->SetMarkerColor(histos_DsToPiKK.back()->GetLineColor());
        histos_DsToPiKK.back()->SetFillColor(histos_DsToPiKK.back()->GetLineColor());
        histos_DsToPiKK.back()->SetFillStyle(3004);
        max_values.push_back(histos_DsToPiKK.back()->GetMaximum());
        //
        histos_DsToPiKKPi0.back()->SetTitle("D_{s}^{+}#rightarrowK^{#minus}#pi^{+}#pi^{+}#pi^{0}");
        histos_DsToPiKKPi0.back()->SetLineColor(colors.at(DsToPiKKPi0-1));
        histos_DsToPiKKPi0.back()->SetLineWidth(2);
        histos_DsToPiKKPi0.back()->SetMarkerColor(histos_DsToPiKKPi0.back()->GetLineColor());
        histos_DsToPiKKPi0.back()->SetFillColor(histos_DsToPiKKPi0.back()->GetLineColor());
        histos_DsToPiKKPi0.back()->SetFillStyle(3005);
        max_values.push_back(histos_DsToPiKKPi0.back()->GetMaximum());
        //
        histos_DsToPiPiK.back()->SetTitle("D_{s}^{+}#rightarrowK^{+}#pi^{+}#pi^{#minus}");
        histos_DsToPiPiK.back()->SetLineColor(colors.at(DsToPiPiK-1));
        histos_DsToPiPiK.back()->SetLineWidth(2);
        histos_DsToPiPiK.back()->SetMarkerColor(histos_DsToPiPiK.back()->GetLineColor());
        histos_DsToPiPiK.back()->SetFillColor(histos_DsToPiPiK.back()->GetLineColor());
        histos_DsToPiPiK.back()->SetFillStyle(3006);
        max_values.push_back(histos_DsToPiPiK.back()->GetMaximum());
        //
        histos_DsToPiPiPi.back()->SetTitle("D_{s}^{+}#rightarrow#pi^{+}#pi^{#minus}#pi^{+}");
        histos_DsToPiPiPi.back()->SetLineColor(colors.at(DsToPiPiPi-1));
        histos_DsToPiPiPi.back()->SetLineWidth(2);
        histos_DsToPiPiPi.back()->SetMarkerColor(histos_DsToPiPiPi.back()->GetLineColor());
        histos_DsToPiPiPi.back()->SetFillColor(histos_DsToPiPiPi.back()->GetLineColor());
        histos_DsToPiPiPi.back()->SetFillStyle(3007);
        max_values.push_back(histos_DsToPiPiPi.back()->GetMaximum());
        //
        histos_DsToPiPiPiPi0.back()->SetTitle("D_{s}^{+}#rightarrow#pi^{+}#pi^{#minus}#pi^{+}#pi^{0}");
        histos_DsToPiPiPiPi0.back()->SetLineColor(colors.at(DsToPiPiPiPi0-1));
        histos_DsToPiPiPiPi0.back()->SetLineWidth(2);
        histos_DsToPiPiPiPi0.back()->SetMarkerColor(histos_DsToPiPiPiPi0.back()->GetLineColor());
        histos_DsToPiPiPiPi0.back()->SetFillColor(histos_DsToPiPiPiPi0.back()->GetLineColor());
        histos_DsToPiPiPiPi0.back()->SetFillStyle(3008);
        max_values.push_back(histos_DsToPiPiPiPi0.back()->GetMaximum());
        // D* channels
        histos_DstarToPiKPi.back()->SetTitle("D^{*+}#rightarrow#pi^{+}K^{#minus}#pi^{+}");
        histos_DstarToPiKPi.back()->SetLineColor(colors.at(DstarToPiKPi-1));
        histos_DstarToPiKPi.back()->SetLineWidth(2);
        histos_DstarToPiKPi.back()->SetMarkerColor(histos_DstarToPiKPi.back()->GetLineColor());
        histos_DstarToPiKPi.back()->SetFillColor(histos_DstarToPiKPi.back()->GetLineColor());
        histos_DstarToPiKPi.back()->SetFillStyle(3004);
        max_values.push_back(histos_DstarToPiKPi.back()->GetMaximum());
        //
        histos_DstarToPiKPiPi0.back()->SetTitle("D^{*+}#rightarrow#pi^{+}K^{#minus}#pi^{+}#pi^{0}");
        histos_DstarToPiKPiPi0.back()->SetLineColor(colors.at(DstarToPiKPiPi0-1));
        histos_DstarToPiKPiPi0.back()->SetLineWidth(2);
        histos_DstarToPiKPiPi0.back()->SetMarkerColor(histos_DstarToPiKPiPi0.back()->GetLineColor());
        histos_DstarToPiKPiPi0.back()->SetFillColor(histos_DstarToPiKPiPi0.back()->GetLineColor());
        histos_DstarToPiKPiPi0.back()->SetFillStyle(3005);
        max_values.push_back(histos_DstarToPiKPiPi0.back()->GetMaximum());
        //
        histos_DstarToPiKPiPi0Pi0.back()->SetTitle("D^{*+}#rightarrow#pi^{+}K^{#minus}#pi^{+}#pi^{0}#pi^{0}");
        histos_DstarToPiKPiPi0Pi0.back()->SetLineColor(colors.at(DstarToPiKPiPi0Pi0-1));
        histos_DstarToPiKPiPi0Pi0.back()->SetLineWidth(2);
        histos_DstarToPiKPiPi0Pi0.back()->SetMarkerColor(histos_DstarToPiKPiPi0Pi0.back()->GetLineColor());
        histos_DstarToPiKPiPi0Pi0.back()->SetFillColor(histos_DstarToPiKPiPi0Pi0.back()->GetLineColor());
        histos_DstarToPiKPiPi0Pi0.back()->SetFillStyle(3006);
        max_values.push_back(histos_DstarToPiKPiPi0Pi0.back()->GetMaximum());
        //
        histos_DstarToPiKK.back()->SetTitle("D^{*+}#rightarrow#pi^{+}K^{#minus}K^{+}");
        histos_DstarToPiKK.back()->SetLineColor(colors.at(DstarToPiKK-1));
        histos_DstarToPiKK.back()->SetLineWidth(2);
        histos_DstarToPiKK.back()->SetMarkerColor(histos_DstarToPiKK.back()->GetLineColor());
        histos_DstarToPiKK.back()->SetFillColor(histos_DstarToPiKK.back()->GetLineColor());
        histos_DstarToPiKK.back()->SetFillStyle(3007);
        max_values.push_back(histos_DstarToPiKK.back()->GetMaximum());
        //
        histos_DstarToPiKKPi0.back()->SetTitle("D^{*+}#rightarrow#pi^{+}K^{#minus}K^{+}#pi^{0}");
        histos_DstarToPiKKPi0.back()->SetLineColor(colors.at(DstarToPiKKPi0-1));
        histos_DstarToPiKKPi0.back()->SetLineWidth(2);
        histos_DstarToPiKKPi0.back()->SetMarkerColor(histos_DstarToPiKKPi0.back()->GetLineColor());
        histos_DstarToPiKKPi0.back()->SetFillColor(histos_DstarToPiKKPi0.back()->GetLineColor());
        histos_DstarToPiKKPi0.back()->SetFillStyle(3008);
        max_values.push_back(histos_DstarToPiKKPi0.back()->GetMaximum());
        //
        histos_DstarToPiPiPi.back()->SetTitle("D^{*+}#rightarrow#pi^{+}#pi^{#minus}#pi^{+}");
        histos_DstarToPiPiPi.back()->SetLineColor(colors.at(DstarToPiPiPi-1));
        histos_DstarToPiPiPi.back()->SetLineWidth(2);
        histos_DstarToPiPiPi.back()->SetMarkerColor(histos_DstarToPiPiPi.back()->GetLineColor());
        histos_DstarToPiPiPi.back()->SetFillColor(histos_DstarToPiPiPi.back()->GetLineColor());
        histos_DstarToPiPiPi.back()->SetFillStyle(3009);
        max_values.push_back(histos_DstarToPiPiPi.back()->GetMaximum());
        //
        histos_DstarToPiPiPiPi0.back()->SetTitle("D^{*+}#rightarrow#pi^{+}#pi^{#minus}#pi^{+}#pi^{0}");
        histos_DstarToPiPiPiPi0.back()->SetLineColor(colors.at(DstarToPiPiPiPi0-1));
        histos_DstarToPiPiPiPi0.back()->SetLineWidth(2);
        histos_DstarToPiPiPiPi0.back()->SetMarkerColor(histos_DstarToPiPiPiPi0.back()->GetLineColor());
        histos_DstarToPiPiPiPi0.back()->SetFillColor(histos_DstarToPiPiPiPi0.back()->GetLineColor());
        histos_DstarToPiPiPiPi0.back()->SetFillStyle(3010);
        max_values.push_back(histos_DstarToPiPiPiPi0.back()->GetMaximum());
        // Lc channels
        histos_LcToPKPi.back()->SetTitle("#Lambda_{c}^{+}#rightarrowpK^{#minus}#pi^{+}");
        histos_LcToPKPi.back()->SetLineColor(colors.at(LcToPKPi-1));
        histos_LcToPKPi.back()->SetLineWidth(2);
        histos_LcToPKPi.back()->SetMarkerColor(histos_LcToPKPi.back()->GetLineColor());
        histos_LcToPKPi.back()->SetFillColor(histos_LcToPKPi.back()->GetLineColor());
        histos_LcToPKPi.back()->SetFillStyle(3001);
        max_values.push_back(histos_LcToPKPi.back()->GetMaximum());
        //
        histos_LcToPKPiPi0.back()->SetTitle("#Lambda_{c}^{+}#rightarrowpK^{#minus}#pi^{+}#pi^{0}");
        histos_LcToPKPiPi0.back()->SetLineColor(colors.at(LcToPKPiPi0-1));
        histos_LcToPKPiPi0.back()->SetLineWidth(2);
        histos_LcToPKPiPi0.back()->SetMarkerColor(histos_LcToPKPiPi0.back()->GetLineColor());
        histos_LcToPKPiPi0.back()->SetFillColor(histos_LcToPKPiPi0.back()->GetLineColor());
        histos_LcToPKPiPi0.back()->SetFillStyle(3004);
        max_values.push_back(histos_LcToPKPiPi0.back()->GetMaximum());
        //
        histos_LcToPPiPi.back()->SetTitle("#Lambda_{c}^{+}#rightarrowp#pi^{#minus}#pi^{+}");
        histos_LcToPPiPi.back()->SetLineColor(colors.at(LcToPPiPi-1));
        histos_LcToPPiPi.back()->SetLineWidth(2);
        histos_LcToPPiPi.back()->SetMarkerColor(histos_LcToPPiPi.back()->GetLineColor());
        histos_LcToPPiPi.back()->SetFillColor(histos_LcToPPiPi.back()->GetLineColor());
        histos_LcToPPiPi.back()->SetFillStyle(3005);
        max_values.push_back(histos_LcToPPiPi.back()->GetMaximum());
        //
        histos_LcToPKK.back()->SetTitle("#Lambda_{c}^{+}#rightarrowpK^{#minus}K^{+}");
        histos_LcToPKK.back()->SetLineColor(colors.at(LcToPKK-1));
        histos_LcToPKK.back()->SetLineWidth(2);
        histos_LcToPKK.back()->SetMarkerColor(histos_LcToPKK.back()->GetLineColor());
        histos_LcToPKK.back()->SetFillColor(histos_LcToPKK.back()->GetLineColor());
        histos_LcToPKK.back()->SetFillStyle(3006);
        max_values.push_back(histos_LcToPKK.back()->GetMaximum());
        // Xic channels
        histos_XicToPKPi.back()->SetTitle("#Xi_{c}^{+}#rightarrowpK^{#minus}#pi^{+}");
        histos_XicToPKPi.back()->SetLineColor(colors.at(XicToPKPi-1));
        histos_XicToPKPi.back()->SetLineWidth(2);
        histos_XicToPKPi.back()->SetMarkerColor(histos_XicToPKPi.back()->GetLineColor());
        histos_XicToPKPi.back()->SetFillColor(histos_XicToPKPi.back()->GetLineColor());
        histos_XicToPKPi.back()->SetFillStyle(3004);
        max_values.push_back(histos_XicToPKPi.back()->GetMaximum());
        //
        histos_XicToPKK.back()->SetTitle("#Xi_{c}^{+}#rightarrowpK^{#minus}K^{+}");
        histos_XicToPKK.back()->SetLineColor(colors.at(XicToPKK-1));
        histos_XicToPKK.back()->SetLineWidth(2);
        histos_XicToPKK.back()->SetMarkerColor(histos_XicToPKK.back()->GetLineColor());
        histos_XicToPKK.back()->SetFillColor(histos_XicToPKK.back()->GetLineColor());
        histos_XicToPKK.back()->SetFillStyle(3005);
        max_values.push_back(histos_XicToPKK.back()->GetMaximum());
        //
        histos_XicToSPiPi.back()->SetTitle("#Xi_{c}^{+}#rightarrow#Sigma^{+}#pi^{#minus}#pi^{+}");
        histos_XicToSPiPi.back()->SetLineColor(colors.at(XicToSPiPi-1));
        histos_XicToSPiPi.back()->SetLineWidth(2);
        histos_XicToSPiPi.back()->SetMarkerColor(histos_XicToSPiPi.back()->GetLineColor());
        histos_XicToSPiPi.back()->SetFillColor(histos_XicToSPiPi.back()->GetLineColor());
        histos_XicToSPiPi.back()->SetFillStyle(3006);
        max_values.push_back(histos_XicToSPiPi.back()->GetMaximum());
        // sum
        histos_bkgSum.back()->SetTitle("sum of backgrounds");
        histos_bkgSum.back()->SetLineColor(kGray+1);
        histos_bkgSum.back()->SetLineWidth(2);
        histos_bkgSum.back()->SetMarkerColor(histos_bkgSum.back()->GetLineColor());

        /// determine the maximum Y
        maxY = *std::max_element(max_values.begin(), max_values.end()) * 1.5;
        maxY = std::max(maxY, 1.f);

        /// Draw stuff
        // D+ channels
        TH1F* frame_DplusChannels = can_DplusChannels->cd(ptBin+1)->DrawFrame(2.1, 0.00000001, 2.5, maxY, title.c_str());
        frame_DplusChannels->GetXaxis()->SetTitle("#it{M}(pK#pi) (GeV/#it{c}^{2})");
        frame_DplusChannels->GetYaxis()->SetTitle(Form("counts/%.0f MeV/#it{c}^{2}", binwidth));
        frame_DplusChannels->SetLineWidth(0);
        histos_DplusToPiKPi.back()->Draw("samehist");
        histos_DplusToPiKPiPi0.back()->Draw("samehist");
        histos_DplusToPiPiPi.back()->Draw("samehist");
        histos_DplusToPiKK.back()->Draw("samehist");
        gPad->Update();
        gPad->SetBottomMargin(0.15);
        gPad->SetLeftMargin(0.15);
        gPad->SetTicks();
        gPad->Update();
        // Ds channels
        TH1F* frame_DsChannels = can_DsChannels->cd(ptBin+1)->DrawFrame(2.1, 0.00000001, 2.5, maxY, title.c_str());
        frame_DsChannels->GetXaxis()->SetTitle("#it{M}(pK#pi) (GeV/#it{c}^{2})");
        frame_DsChannels->GetYaxis()->SetTitle(Form("counts/%.0f MeV/#it{c}^{2}", binwidth));
        frame_DsChannels->SetLineWidth(0);
        histos_DsToPiKK.back()->Draw("samehist");
        histos_DsToPiKKPi0.back()->Draw("samehist");
        histos_DsToPiPiK.back()->Draw("samehist");
        histos_DsToPiPiPi.back()->Draw("samehist");
        histos_DsToPiPiPiPi0.back()->Draw("samehist");
        gPad->Update();
        gPad->SetBottomMargin(0.15);
        gPad->SetLeftMargin(0.15);
        gPad->SetTicks();
        gPad->Update();
        // D* channels
        TH1F* frame_DstarChannels = can_DstarChannels->cd(ptBin+1)->DrawFrame(2.1, 0.00000001, 2.5, maxY, title.c_str());
        frame_DstarChannels->GetXaxis()->SetTitle("#it{M}(pK#pi) (GeV/#it{c}^{2})");
        frame_DstarChannels->GetYaxis()->SetTitle(Form("counts/%.0f MeV/#it{c}^{2}", binwidth));
        frame_DstarChannels->SetLineWidth(0);
        histos_DstarToPiKPi.back()->Draw("samehist");
        histos_DstarToPiKPiPi0.back()->Draw("samehist");
        histos_DstarToPiKPiPi0Pi0.back()->Draw("samehist");
        histos_DstarToPiKK.back()->Draw("samehist");
        histos_DstarToPiKKPi0.back()->Draw("samehist");
        histos_DstarToPiPiPi.back()->Draw("samehist");
        histos_DstarToPiPiPiPi0.back()->Draw("samehist");
        gPad->Update();
        gPad->SetBottomMargin(0.15);
        gPad->SetLeftMargin(0.15);
        gPad->SetTicks();
        gPad->Update();
        // Lc channels
        TH1F* frame_LcChannels = can_LcChannels->cd(ptBin+1)->DrawFrame(2.1, 0.00000001, 2.5, maxY, title.c_str());
        frame_LcChannels->GetXaxis()->SetTitle("#it{M}(pK#pi) (GeV/#it{c}^{2})");
        frame_LcChannels->GetYaxis()->SetTitle(Form("counts/%.0f MeV/#it{c}^{2}", binwidth));
        frame_LcChannels->SetLineWidth(0);
        histos_LcToPKPi.back()->Draw("samehist");
        histos_LcToPKPiPi0.back()->Draw("samehist");
        histos_LcToPPiPi.back()->Draw("samehist");
        histos_LcToPKK.back()->Draw("samehist");
        gPad->Update();
        gPad->SetBottomMargin(0.15);
        gPad->SetLeftMargin(0.15);
        gPad->SetTicks();
        gPad->Update();
        // Xic channel
        TH1F* frame_XicChannels = can_XicChannels->cd(ptBin+1)->DrawFrame(2.1, 0.00000001, 2.5, maxY, title.c_str());
        frame_XicChannels->GetXaxis()->SetTitle("#it{M}(pK#pi) (GeV/#it{c}^{2})");
        frame_XicChannels->GetYaxis()->SetTitle(Form("counts/%.0f MeV/#it{c}^{2}", binwidth));
        frame_XicChannels->SetLineWidth(0);
        histos_XicToPKPi.back()->Draw("samehist");
        histos_XicToPKK.back()->Draw("samehist");
        histos_XicToSPiPi.back()->Draw("samehist");
        gPad->Update();
        gPad->SetBottomMargin(0.15);
        gPad->SetLeftMargin(0.15);
        gPad->SetTicks();
        gPad->Update();
        // all channels
        TH1F* frame_allChannels = can_allChannels->cd(ptBin+1)->DrawFrame(2.1, 0.00000001, 2.5, maxY, title.c_str());
        frame_allChannels->GetXaxis()->SetTitle("#it{M}(pK#pi) (GeV/#it{c}^{2})");
        frame_allChannels->GetYaxis()->SetTitle(Form("counts/%.0f MeV/#it{c}^{2}", binwidth));
        histos_DplusToPiKPi.back()->Draw("samehist");
        histos_DplusToPiKPiPi0.back()->Draw("samehist");
        histos_DplusToPiPiPi.back()->Draw("samehist");
        histos_DplusToPiKK.back()->Draw("samehist");
        histos_DsToPiKK.back()->Draw("samehist");
        histos_DsToPiKKPi0.back()->Draw("samehist");
        histos_DsToPiPiK.back()->Draw("samehist");
        histos_DsToPiPiPi.back()->Draw("samehist");
        histos_DsToPiPiPiPi0.back()->Draw("samehist");
        histos_DstarToPiKPi.back()->Draw("samehist");
        histos_DstarToPiKPiPi0.back()->Draw("samehist");
        histos_DstarToPiKPiPi0Pi0.back()->Draw("samehist");
        histos_DstarToPiKK.back()->Draw("samehist");
        histos_DstarToPiKKPi0.back()->Draw("samehist");
        histos_DstarToPiPiPi.back()->Draw("samehist");
        histos_DstarToPiPiPiPi0.back()->Draw("samehist");
        histos_LcToPKPi.back()->Draw("samehist");
        histos_LcToPKPiPi0.back()->Draw("samehist");
        histos_LcToPPiPi.back()->Draw("samehist");
        histos_LcToPKK.back()->Draw("samehist");
        histos_XicToPKPi.back()->Draw("samehist");
        histos_XicToPKK.back()->Draw("samehist");
        histos_XicToSPiPi.back()->Draw("samehist");
        histos_bkgSum.back()->Draw("samehist");
        gPad->Update();
        gPad->SetBottomMargin(0.15);
        gPad->SetLeftMargin(0.15);
        gPad->SetTicks();
        gPad->Update();

        
        if(ptBin == row*col-2 || ptBin == row*col-1) {
            int ptBin_legend = ptBin;
            float minY_legend = 0.6;
            float maxY_legend = 0.8;
            if (ptBin == row*col-2) {
                ptBin_legend++;
                minY_legend = 0.2;
                maxY_legend = 0.85;
            }
            // D+ channels
            TLegend* leg_DplusChannels = new TLegend(0.2, minY_legend, 0.8, maxY_legend);
            leg_DplusChannels->SetFillStyle(0);
            leg_DplusChannels->SetBorderSize(0);
            leg_DplusChannels->AddEntry(histos_DplusToPiKPi.back());
            leg_DplusChannels->AddEntry(histos_DplusToPiKPiPi0.back());
            leg_DplusChannels->AddEntry(histos_DplusToPiPiPi.back());
            leg_DplusChannels->AddEntry(histos_DplusToPiKK.back());
            can_DplusChannels->cd(ptBin_legend+1);
            leg_DplusChannels->Draw();
            gPad->Update();
            // Ds channels
            TLegend* leg_DsChannels = new TLegend(0.2, minY_legend, 0.8, maxY_legend);
            leg_DsChannels->SetFillStyle(0);
            leg_DsChannels->SetBorderSize(0);
            leg_DsChannels->AddEntry(histos_DsToPiKK.back());
            leg_DsChannels->AddEntry(histos_DsToPiKKPi0.back());
            leg_DsChannels->AddEntry(histos_DsToPiPiK.back());
            leg_DsChannels->AddEntry(histos_DsToPiPiPi.back());
            leg_DsChannels->AddEntry(histos_DsToPiPiPiPi0.back());
            can_DsChannels->cd(ptBin_legend+1);
            leg_DsChannels->Draw();
            gPad->Update();
            // D* channels
            TLegend* leg_DstarChannels = new TLegend(0.2, minY_legend, 0.8, maxY_legend);
            leg_DstarChannels->SetFillStyle(0);
            leg_DstarChannels->SetBorderSize(0);
            leg_DstarChannels->AddEntry(histos_DstarToPiKPi.back());
            leg_DstarChannels->AddEntry(histos_DstarToPiKPiPi0.back());
            leg_DstarChannels->AddEntry(histos_DstarToPiKPiPi0Pi0.back());
            leg_DstarChannels->AddEntry(histos_DstarToPiKK.back());
            leg_DstarChannels->AddEntry(histos_DstarToPiKKPi0.back());
            leg_DstarChannels->AddEntry(histos_DstarToPiPiPi.back());
            leg_DstarChannels->AddEntry(histos_DstarToPiPiPiPi0.back());
            can_DstarChannels->cd(ptBin_legend+1);
            leg_DstarChannels->Draw();
            gPad->Update();
            // Lc channels
            TLegend* leg_LcChannels = new TLegend(0.2, minY_legend, 0.8, maxY_legend);
            leg_LcChannels->SetFillStyle(0);
            leg_LcChannels->SetBorderSize(0);
            leg_LcChannels->AddEntry(histos_LcToPKPi.back());
            leg_LcChannels->AddEntry(histos_LcToPKPiPi0.back());
            leg_LcChannels->AddEntry(histos_LcToPPiPi.back());
            leg_LcChannels->AddEntry(histos_LcToPKK.back());
            can_LcChannels->cd(ptBin_legend+1);
            leg_LcChannels->Draw();
            gPad->Update();
            // Xic channels
            TLegend* leg_XicChannels = new TLegend(0.2, minY_legend, 0.8, maxY_legend);
            leg_XicChannels->SetFillStyle(0);
            leg_XicChannels->SetBorderSize(0);
            leg_XicChannels->AddEntry(histos_XicToPKPi.back());
            leg_XicChannels->AddEntry(histos_XicToPKK.back());
            leg_XicChannels->AddEntry(histos_XicToSPiPi.back());
            can_XicChannels->cd(ptBin_legend+1);
            leg_XicChannels->Draw();
            gPad->Update();
            // all channels
            TLegend* leg_allChannels = new TLegend(0.2, minY_legend, 0.8, maxY_legend);
            leg_allChannels->SetFillStyle(0);
            leg_allChannels->SetBorderSize(0);
            leg_allChannels->SetNColumns(2);
            leg_allChannels->AddEntry(histos_DplusToPiKPi.back());
            leg_allChannels->AddEntry(histos_DplusToPiKPiPi0.back());
            leg_allChannels->AddEntry(histos_DplusToPiPiPi.back());
            leg_allChannels->AddEntry(histos_DplusToPiKK.back());
            leg_allChannels->AddEntry(histos_DsToPiKK.back());
            leg_allChannels->AddEntry(histos_DsToPiKKPi0.back());
            leg_allChannels->AddEntry(histos_DsToPiPiK.back());
            leg_allChannels->AddEntry(histos_DsToPiPiPi.back());
            leg_allChannels->AddEntry(histos_DsToPiPiPiPi0.back());
            leg_allChannels->AddEntry(histos_DstarToPiKPi.back());
            leg_allChannels->AddEntry(histos_DstarToPiKPiPi0.back());
            leg_allChannels->AddEntry(histos_DstarToPiKPiPi0Pi0.back());
            leg_allChannels->AddEntry(histos_DstarToPiKK.back());
            leg_allChannels->AddEntry(histos_DstarToPiKKPi0.back());
            leg_allChannels->AddEntry(histos_DstarToPiPiPi.back());
            leg_allChannels->AddEntry(histos_DstarToPiPiPiPi0.back());
            leg_allChannels->AddEntry(histos_LcToPKPi.back());
            leg_allChannels->AddEntry(histos_LcToPKPiPi0.back());
            leg_allChannels->AddEntry(histos_LcToPPiPi.back());
            leg_allChannels->AddEntry(histos_LcToPKK.back());
            leg_allChannels->AddEntry(histos_XicToPKPi.back());
            leg_allChannels->AddEntry(histos_XicToPKK.back());
            leg_allChannels->AddEntry(histos_XicToSPiPi.back());
            leg_allChannels->AddEntry(histos_bkgSum.back());
            can_allChannels->cd(ptBin_legend+1);
            leg_allChannels->Draw();
            gPad->Update();
        }
        

    } /// end loop over pt intervals

    /// save stuff
    TFile* fout = new TFile(Form("file_corrBkgPKPi_MC%s%s.root", scaleByBrs ? "_brScaled" : "", applyBdt ? "_bdtCuts" : ""), "recreate");
    can_allChannels->Write();
    for(int ptBin=0; ptBin<ptMins.size(); ptBin++) {
        fout ->cd();
        TDirectory* td = fout->mkdir(Form("ptBin%d", ptBin));
        td->cd();
        histos_DplusToPiKPi.at(ptBin)->Write();
        histos_DplusToPiKPiPi0.at(ptBin)->Write();
        histos_DplusToPiPiPi.at(ptBin)->Write();
        histos_DplusToPiKK.at(ptBin)->Write();
        histos_DsToPiKK.at(ptBin)->Write();
        histos_DsToPiKKPi0.at(ptBin)->Write();
        histos_DsToPiPiK.at(ptBin)->Write();
        histos_DsToPiPiPi.at(ptBin)->Write();
        histos_DsToPiPiPiPi0.at(ptBin)->Write();
        histos_DstarToPiKPi.at(ptBin)->Write();
        histos_DstarToPiKPiPi0.at(ptBin)->Write();
        histos_DstarToPiKPiPi0Pi0.at(ptBin)->Write();
        histos_DstarToPiKK.at(ptBin)->Write();
        histos_DstarToPiKKPi0.at(ptBin)->Write();
        histos_DstarToPiPiPi.at(ptBin)->Write();
        histos_DstarToPiPiPiPi0.at(ptBin)->Write();
        histos_LcToPKPi.at(ptBin)->Write();
        histos_LcToPKPiPi0.at(ptBin)->Write();
        histos_LcToPPiPi.at(ptBin)->Write();
        histos_LcToPKK.at(ptBin)->Write();
        histos_XicToPKPi.at(ptBin)->Write();
        histos_XicToPKK.at(ptBin)->Write();
        histos_XicToSPiPi.at(ptBin)->Write();
        histos_bkgSum.at(ptBin)->Write();
    }

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
