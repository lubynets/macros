//
// Created by oleksii on 12.09.2025.
//

#include "HelperGeneral.hpp"
#include "HelperMath.hpp"

#include <TAxis.h>
#include <THnSparse.h>
#include <TH1.h>
#include <TH2.h>

#include <iostream>
#include <string>
#include <string_view>
#include <vector>

using namespace HelperGeneral;
using namespace HelperMath;
using namespace std::string_literals;

bool gIsDoWeight{false};
std::vector<double> gBdtSignalLowerValues{};

std::vector<double> pTRanges = {3, 4, 5, 8, 12, 20};
const std::vector<double> bdtBgUpperValuesVsPt = {0.02, 0.02, 0.02, 0.04, 0.08}; // standard
// const std::vector<double> bdtBgUpperValuesVsPt = {0.01, 0.01, 0.01, 0.02, 0.04}; // tight
// const std::vector<double> bdtBgUpperValuesVsPt = {0.03, 0.03, 0.03, 0.06, 0.12}; // loose
// const std::vector<double> bdtBgUpperValuesVsPt = {0.015, 0.015, 0.015, 0.03, 0.06}; // semitight
// const std::vector<double> bdtBgUpperValuesVsPt = {0.05, 0.05, 0.05, 0.08, 0.16}; // veryloose
// const std::vector<double> bdtBgUpperValuesVsPt = {0.025, 0.025, 0.025, 0.05, 0.10}; // semiloose

constexpr std::string_view LifetimeAxisTitle = /*"T_{proper} (ps)"*/"#it{t}_{proper} (ps)";
constexpr std::string_view LifetimeGenAxisTitle = "#it{t}_{proper, gen} (ps)";
constexpr std::string_view PtAxisTitle = "#it{p}_{T}(#Lambda_{c}^{+}) (GeV/#it{c})";
constexpr std::string_view PtBAxisTitle = "#it{p}_{T}^{B} (GeV/#it{c})";
constexpr std::string_view BgAxisTitle = "BDT bkg score (Lc)";
constexpr std::string_view NpAxisTitle = "BDT non-prompt score (Lc)";
constexpr std::string_view SignalTypeAxisTitle = "candidates type";

const std::string fileOutName{"yield_lifetime_qa_thn.root"};
constexpr bool IsVerbose{true};

const double tauPythia{0.2005};
const double tauPdg{0.2026};
const double sigmaPdg{0.001};

const std::vector<std::pair<std::string, double>> promptnesses {
  {"prompt", 1.},
  {"nonprompt", 2.}
};
enum : int {
  RunOnly = 0,
  RunAndMerge,
  MergeOnly,
  NModeRuns
};
const std::array<std::string, 2> weightsPresences{"", "_W"};

std::string GetPtCutName(size_t iPt);

void FillYield(const std::string& fileName, const std::string& filePtWeightName, const bool isRec, const double nCtSigmaFromPdg, const std::string& pCMM, const std::string& npCMM) {
  const double tauTarget = tauPdg + nCtSigmaFromPdg*sigmaPdg;

  TF1* ctWeight = new TF1("ctWeight", "std::exp(-[0]*x)", 0., 2.);
  ctWeight->SetParameter(0, (tauPythia - tauTarget)/tauTarget/tauPythia);

  if(IsVerbose) std::cout << "FillYield() started\n";
  TFile* fileIn = OpenFileWithNullptrCheck(fileName);
  const std::string fileOpenOption = isRec ? "recreate" : "update";
  TFile* fileOut = TFile::Open(fileOutName.c_str(), fileOpenOption.c_str());
  TFile* fileWeight = gIsDoWeight ? OpenFileWithNullptrCheck(filePtWeightName) : nullptr;
  TH1* histoWeightPrompt = gIsDoWeight ? GetObjectWithNullptrCheck<TH1>(fileWeight, "histoWeight" + pCMM) : nullptr;
  TH1* histoWeightNonPrompt = gIsDoWeight ? GetObjectWithNullptrCheck<TH1>(fileWeight, "histoNPWeight" + npCMM) : nullptr;
  THnSparse* histoRecOrGen = GetObjectWithNullptrCheck<THnSparse>(fileIn, "hf-task-lc/"s + (isRec ? "hnLcVarsWithBdt" : "hnLcVarsGen"));
  const std::map<std::string_view, int> axesIndices = MapTHnSparseAxesIndices(histoRecOrGen);
  THnSparse* histoPromptWeighted = gIsDoWeight ? dynamic_cast<THnSparse*>(histoRecOrGen->Clone()) : nullptr;
  THnSparse* histoNonPromptWeighted = gIsDoWeight ? dynamic_cast<THnSparse*>(histoRecOrGen->Clone()) : nullptr;
  if (gIsDoWeight) {
    ScaleTHnSparseWithWeight(histoPromptWeighted, axesIndices.at(PtAxisTitle), histoWeightPrompt);
    ScaleTHnSparseWithWeight(histoNonPromptWeighted, axesIndices.at(PtBAxisTitle), histoWeightNonPrompt);

    if(nCtSigmaFromPdg != UndefValueDouble) ScaleTHnSparseWithWeight(histoPromptWeighted, axesIndices.at(LifetimeAxisTitle), ctWeight);
  }
  
  auto ProcessTHnSparse = [&](THnSparse* histoIn, const std::string& histoNameSuffix="", const std::vector<std::pair<std::string, double>>& promptnessesToProcess=promptnesses) {
    if(IsVerbose) std::cout << "ProcessTHnSparse() started\n";
    CheckTAxisForRanges(*histoIn->GetAxis(axesIndices.at(PtAxisTitle)), pTRanges);
    CheckTAxisForRanges(*histoIn->GetAxis(axesIndices.at(SignalTypeAxisTitle)), {1., 2., 3.});
    if(isRec) CheckTAxisForRanges(*histoIn->GetAxis(axesIndices.at(BgAxisTitle)), bdtBgUpperValuesVsPt);

    for(size_t iPt=0, nPts=pTRanges.size()-1; iPt<nPts; ++iPt) {
      if(IsVerbose) std::cout << "ProcessTHnSparse(): iPt = " << iPt << "\n";
      SetTHnSparseAxisRanges(histoIn, axesIndices.at(PtAxisTitle), pTRanges.at(iPt), pTRanges.at(iPt + 1));
      if(isRec) SetTHnSparseAxisRanges(histoIn, axesIndices.at(BgAxisTitle), 0., bdtBgUpperValuesVsPt.at(iPt));
      for(const auto& promptness : promptnessesToProcess) {
        if(IsVerbose) std::cout << "ProcessTHnSparse(): promptness = " << promptness.first << "\n";
        const std::string dirName = (isRec ? "rec/" : "gen/") + promptness.first + "/" + GetPtCutName(iPt);
        SetTHnSparseAxisRanges(histoIn, axesIndices.at(SignalTypeAxisTitle), promptness.second, promptness.second+1.f);
        // for rec - real gBdtSignalLowerValues; for gen - fake 1-element vector for universality reasons
        const auto& bdtSignalLowerValues = isRec ? gBdtSignalLowerValues : std::vector<double>{UndefValueDouble};
        if(IsVerbose) std::cout << "ProcessTHnSparse(): bsc = ";
        for(const auto& bsc : bdtSignalLowerValues) {
          if(IsVerbose) std::cout << bsc << "\t";
          const std::string histoName = isRec ?
                                        "hT_NPgt" + to_string_with_precision(bsc, 2) + histoNameSuffix :
                                        "hT" + histoNameSuffix;

          if(isRec) SetTHnSparseAxisRanges(histoIn, axesIndices.at(NpAxisTitle), bsc, 1.);
          TH1* histoYield = histoIn->Projection(axesIndices.at(LifetimeAxisTitle));
          histoYield->SetDirectory(nullptr);
          CD(fileOut, dirName);
          histoYield->Write(histoName.c_str());
          if(isRec) {
            std::string histoName2D{histoName};
            ReplaceSubstrInStr(histoName2D, "hT", "h2T");
            TH2* histoYield2D = histoIn->Projection(axesIndices.at(LifetimeGenAxisTitle), axesIndices.at(LifetimeAxisTitle));
            histoYield2D->SetDirectory(nullptr);
            CD(fileOut, dirName);
            histoYield2D->Write(histoName2D.c_str());
          }
          if(isRec) SetTHnSparseAxisRanges(histoIn, axesIndices.at(NpAxisTitle));
        } // bdtSignalLowerValues
        if(IsVerbose) std::cout << "\n";
        SetTHnSparseAxisRanges(histoIn, axesIndices.at(SignalTypeAxisTitle));
      } // promptnessesToProcess
      SetTHnSparseAxisRanges(histoIn, axesIndices.at(PtAxisTitle));
      if(isRec) SetTHnSparseAxisRanges(histoIn, axesIndices.at(BgAxisTitle));
    } // pTRanges
    if(IsVerbose) std::cout << "ProcessTHnSparse() finished\n";
  };
  
  ProcessTHnSparse(histoRecOrGen);
  if (gIsDoWeight) {
    ProcessTHnSparse(histoPromptWeighted, "_W", {promptnesses.at(0)});
    ProcessTHnSparse(histoNonPromptWeighted, "_W", {promptnesses.at(1)});
  }

  fileOut->Close();
  fileIn->Close();
  if (gIsDoWeight) fileWeight->Close();
  if(IsVerbose) std::cout << "FillYield() finished\n";
}

std::string GetPtCutName(size_t iPt) {
  std::pair<size_t, size_t> iPTMinMax = (iPt == pTRanges.size()-1) ? std::pair<size_t, size_t>{0, pTRanges.size()-1} : std::pair<size_t, size_t>{iPt, iPt+1};
  return "pT_" + to_string_with_precision(pTRanges.at(iPTMinMax.first), 0) + "_" + to_string_with_precision(pTRanges.at(iPTMinMax.second), 0);
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./yield_lifetime_qa_thn fileNameIn (modeRun=RunOnly=0 [RunAndMerge=1, MergeOnly=2]) (filePtWeightName) (nCtSigmaFromPdg)" << std::endl;
    exit(EXIT_FAILURE);
  }
  if(bdtBgUpperValuesVsPt.size() != pTRanges.size() - 1) throw std::runtime_error("bdtUpperValuesVsPt.size() != pTRanges.size() - 1");

  const std::string fileNameIn = argv[1];
  const int modeRun = argc > 2 ? std::stoi(argv[2]) : RunOnly;
  const std::string filePtWeightName = argc > 3 ? argv[3] : "";
  const std::string pCMM = argc> 4 ? argv[4] : "cent";
  const std::string npCMM = argc > 5 ? argv[5] : "cent";

  const double nCtSigmaFromPdg = UndefValueDouble;

  if(modeRun < 0 || modeRun >= NModeRuns) throw std::runtime_error("modeRun < 0 || modeRun >= NModeRuns");

  gIsDoWeight = !filePtWeightName.empty();

  const std::string& fileName = modeRun != MergeOnly ? ReadNthLine(fileNameIn) : fileNameIn;

//   for (int iB = 0; iB <= 99; iB++) {
//     gBdtSignalLowerValues.emplace_back(0.01 * iB);
//   }
  gBdtSignalLowerValues = {0.20, 0.25/*, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90*/};

  if(modeRun != MergeOnly) {
    FillYield(fileName, filePtWeightName, true, nCtSigmaFromPdg, pCMM, npCMM);
    FillYield(fileName, filePtWeightName, false, nCtSigmaFromPdg, pCMM, npCMM);
  }

  if(modeRun == RunOnly) return 0;

  std::vector<std::string> pTCutNames;
  for(size_t iPt=0, nPts=pTRanges.size()-1; iPt<nPts; ++iPt) {
    pTCutNames.emplace_back(GetPtCutName(iPt));
  }

  const double mergePtFrom{3.};
  const double mergePtTo{20.};

  const int skipFirstNBins = std::distance(pTRanges.begin(), std::find(pTRanges.begin(), pTRanges.end(), mergePtFrom));
  const int skipLastNBins = std::distance(std::find(pTRanges.begin(), pTRanges.end(), mergePtTo), pTRanges.end()) - 1;

  if(skipFirstNBins == pTRanges.size() || skipLastNBins == -1) throw std::runtime_error("main(): mergePtFrom or mergePtTo does not correspond to any of pTRanges");

  pTRanges.erase(pTRanges.end()-skipLastNBins, pTRanges.end());
  pTCutNames.erase(pTCutNames.end()-skipLastNBins, pTCutNames.end());
  pTRanges.erase(pTRanges.begin(), pTRanges.begin()+skipFirstNBins);
  pTCutNames.erase(pTCutNames.begin(), pTCutNames.begin()+skipFirstNBins);

  const std::string& mergedFileOutName = modeRun != MergeOnly ? fileOutName : fileName;
  TFile* fileOut = OpenFileWithNullptrCheck(mergedFileOutName.c_str(), "update");

  auto ProcessMerge = [&](const bool isRec) {
    for(const auto& promptness : promptnesses) {
      for(const auto& weightPresence : weightsPresences) {
        if(!gIsDoWeight && weightPresence == "_W") continue;
        // for rec - real gBdtSignalLowerValues; for gen - fake 1-element vector for universality reasons
        const auto& bdtSignalLowerValues = isRec ? gBdtSignalLowerValues : std::vector<double>{UndefValueDouble};
        for (const auto& bslv : bdtSignalLowerValues) {
          std::vector<std::string> histoNames;
          histoNames.reserve(pTCutNames.size());
          for (const auto& ptcn : pTCutNames) {
            histoNames.emplace_back(isRec ?
                                    "rec/" + promptness.first + "/" + ptcn + "/hT_NPgt" + to_string_with_precision(bslv, 2) + weightPresence :
                                    "gen/" + promptness.first + "/" + ptcn + "/hT" + weightPresence);
          } // pTCutNames
          TH1* histoMerged = MergeHistograms<TH1>(fileOut, histoNames);
          CD(fileOut, (isRec ? "rec/" : "gen/") + promptness.first + "/" + GetPtCutName(pTRanges.size() - 1));
          histoMerged->Write((isRec ? "hT_NPgt" + to_string_with_precision(bslv, 2) + weightPresence : "hT" + weightPresence).c_str());
          if(isRec) {
            for(auto& hN : histoNames) {
              ReplaceSubstrInStr(hN, "hT", "h2T");
            }
            TH2* histoMerged = MergeHistograms<TH2>(fileOut, histoNames);
            CD(fileOut, "rec/" + promptness.first + "/" + GetPtCutName(pTRanges.size() - 1));
            histoMerged->Write(("h2T_NPgt" + to_string_with_precision(bslv, 2) + weightPresence).c_str());
          }
        } // bdtSignalLowerValues
      } // weightPresences
    } // promptnesses
  };

  ProcessMerge(true);
  ProcessMerge(false);

  fileOut->Close();

  return 0;
}
