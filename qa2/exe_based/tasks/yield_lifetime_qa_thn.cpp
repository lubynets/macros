//
// Created by oleksii on 12.09.2025.
//

#include "HelperGeneral.hpp"
#include "HelperMath.hpp"

#include <TAxis.h>
#include <THnSparse.h>
#include <TH1.h>

#include <iostream>
#include <string>
#include <string_view>
#include <vector>

using namespace HelperGeneral;
using namespace HelperMath;
using namespace std::string_literals;

bool gIsDoWeight{false};
std::vector<double> gBdtSignalLowerValues{};

std::vector<double> pTRanges = {1, 2, 3, 4, 5, 8, 12, 20};
const std::vector<double> bdtBgUpperValuesVsPt = {0.02, 0.02, 0.02, 0.02, 0.02, 0.04, 0.08};

const std::string_view lifetimeAxisTitle = "T_{proper} (ps)";
const std::string_view pTAxisTitle = "#it{p}_{T}(#Lambda_{c}^{+}) (GeV/#it{c})";
const std::string_view pTBAxisTitle = "#it{p}_{T}^{B} (GeV/#it{c})";
const std::string_view bgAxisTitle = "BDT bkg score (Lc)";
const std::string_view npAxisTitle = "BDT non-prompt score (Lc)";
const std::string_view signalTypeAxisTitle = "candidates type";

const std::string fileOutName{"yield_lifetime_qa_thn.root"};
constexpr bool IsVerbose{true};

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

void FillYield(const std::string& fileName, const std::string& filePtWeightName, const bool isRec) {
  if(IsVerbose) std::cout << "FillYield() started\n";
  TFile* fileIn = OpenFileWithNullptrCheck(fileName);
  const std::string fileOpenOption = isRec ? "recreate" : "update";
  TFile* fileOut = TFile::Open(fileOutName.c_str(), fileOpenOption.c_str());
  TFile* fileWeight = gIsDoWeight ? OpenFileWithNullptrCheck(filePtWeightName) : nullptr;
  TH1* histoWeightPrompt = gIsDoWeight ? GetObjectWithNullptrCheck<TH1>(fileWeight, "histoWeight_pT_0_20") : nullptr;
  TH1* histoWeightNonPrompt = gIsDoWeight ? GetObjectWithNullptrCheck<TH1>(fileWeight, "histoNPWeight") : nullptr;
  THnSparse* histoRecOrGen = GetObjectWithNullptrCheck<THnSparse>(fileIn, "hf-task-lc/"s + (isRec ? "hnLcVarsWithBdt" : "hnLcVarsGen"));
  const std::map<std::string_view, int> axesIndices = MapTHnSparseAxesIndices(histoRecOrGen);
  THnSparse* histoPromptWeighted = gIsDoWeight ? dynamic_cast<THnSparse*>(histoRecOrGen->Clone()) : nullptr;
  THnSparse* histoNonPromptWeighted = gIsDoWeight ? dynamic_cast<THnSparse*>(histoRecOrGen->Clone()) : nullptr;
  if (gIsDoWeight) {
    ScaleTHnSparseWithWeight(histoPromptWeighted, axesIndices.at(pTAxisTitle), histoWeightPrompt);
    ScaleTHnSparseWithWeight(histoNonPromptWeighted, axesIndices.at(pTBAxisTitle), histoWeightNonPrompt);
  }
  
  auto ProcessTHnSparse = [&](THnSparse* histoIn, const std::string& histoNameSuffix="", const std::vector<std::pair<std::string, double>>& promptnessesToProcess=promptnesses) {
    if(IsVerbose) std::cout << "ProcessTHnSparse() started\n";
    CheckTAxisForRanges(*histoIn->GetAxis(axesIndices.at(pTAxisTitle)), pTRanges);
    CheckTAxisForRanges(*histoIn->GetAxis(axesIndices.at(signalTypeAxisTitle)), {1., 2., 3.});
    if(isRec) CheckTAxisForRanges(*histoIn->GetAxis(axesIndices.at(bgAxisTitle)), bdtBgUpperValuesVsPt);

    for(size_t iPt=0, nPts=pTRanges.size()-1; iPt<nPts; ++iPt) {
      if(IsVerbose) std::cout << "ProcessTHnSparse(): iPt = " << iPt << "\n";
      SetTHnSparseAxisRanges(histoIn, axesIndices.at(pTAxisTitle), pTRanges.at(iPt), pTRanges.at(iPt + 1));
      if(isRec) SetTHnSparseAxisRanges(histoIn, axesIndices.at(bgAxisTitle), 0., bdtBgUpperValuesVsPt.at(iPt));
      for(const auto& promptness : promptnessesToProcess) {
        if(IsVerbose) std::cout << "ProcessTHnSparse(): promptness = " << promptness.first << "\n";
        const std::string dirName = (isRec ? "rec/" : "gen/") + promptness.first + "/" + GetPtCutName(iPt);
        SetTHnSparseAxisRanges(histoIn, axesIndices.at(signalTypeAxisTitle), promptness.second, promptness.second+1.f);
        // for rec - real gBdtSignalLowerValues; for gen - fake 1-element vector for universality reasons
        const auto& bdtSignalLowerValues = isRec ? gBdtSignalLowerValues : std::vector<double>{UndefValueDouble};
        if(IsVerbose) std::cout << "ProcessTHnSparse(): bsc = ";
        for(const auto& bsc : bdtSignalLowerValues) {
          if(IsVerbose) std::cout << bsc << "\t";
          const std::string histoName = isRec ?
                                        "hT_NPgt" + to_string_with_precision(bsc, 2) + histoNameSuffix :
                                        "hT" + histoNameSuffix;
          if(isRec) SetTHnSparseAxisRanges(histoIn, axesIndices.at(npAxisTitle), bsc, 1.);
          TH1* histoYield = histoIn->Projection(axesIndices.at(lifetimeAxisTitle));
          histoYield->SetDirectory(nullptr);
          CD(fileOut, dirName);
          histoYield->Write(histoName.c_str());
          if(isRec) SetTHnSparseAxisRanges(histoIn, axesIndices.at(npAxisTitle));
        } // bdtSignalLowerValues
        if(IsVerbose) std::cout << "\n";
        SetTHnSparseAxisRanges(histoIn, axesIndices.at(signalTypeAxisTitle));
      } // promptnessesToProcess
      SetTHnSparseAxisRanges(histoIn, axesIndices.at(pTAxisTitle));
      if(isRec) SetTHnSparseAxisRanges(histoIn, axesIndices.at(bgAxisTitle));
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
    std::cout << " ./yield_lifetime_qa_thn fileNameIn (modeRun=RunOnly=0 [RunAndMerge=1, MergeOnly=2]) (filePtWeightName)" << std::endl;
    exit(EXIT_FAILURE);
  }
  if(bdtBgUpperValuesVsPt.size() != pTRanges.size() - 1) throw std::runtime_error("bdtUpperValuesVsPt.size() != pTRanges.size() - 1");

  const std::string fileNameIn = argv[1];
  const int modeRun = argc > 2 ? std::stoi(argv[2]) : RunOnly;
  const std::string filePtWeightName = argc > 3 ? argv[3] : "";

  if(modeRun < 0 || modeRun >= NModeRuns) throw std::runtime_error("modeRun < 0 || modeRun >= NModeRuns");

  gIsDoWeight = !filePtWeightName.empty();

  const std::string& fileName = modeRun != MergeOnly ? ReadNthLine(fileNameIn) : fileNameIn;

  for (int iB = 0; iB <= 99; iB++) {
    gBdtSignalLowerValues.emplace_back(0.01 * iB);
  }

  if(modeRun != MergeOnly) {
    FillYield(fileName, filePtWeightName, true);
    FillYield(fileName, filePtWeightName, false);
  }

  if(modeRun == RunOnly) return 0;

  const int nLowerPtBinsToExclude{3};
  pTRanges.erase(pTRanges.begin(), pTRanges.begin()+nLowerPtBinsToExclude);

  std::vector<std::string> pTCutNames;
  for(size_t iPt=0, nPts=pTRanges.size()-1; iPt<nPts; ++iPt) {
    pTCutNames.emplace_back(GetPtCutName(iPt));
  }

  const std::string& mergedFileOutName = modeRun != MergeOnly ? fileOutName : fileName;
  TFile* fileOut = OpenFileWithNullptrCheck(mergedFileOutName.c_str(), "update");

  auto ProcessMerge = [&](const bool isRec) {
    for(const auto& promptness : promptnesses) {
      for(const auto& weightPresence : weightsPresences) {
        if(!gIsDoWeight) continue;
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
          TH1* histoMerged = MergeHistograms(fileOut, histoNames);
          CD(fileOut, (isRec ? "rec/" : "gen/") + promptness.first + "/" + GetPtCutName(pTRanges.size() - 1));
          histoMerged->Write((isRec ? "hT_NPgt" + to_string_with_precision(bslv, 2) + weightPresence : "hT" + weightPresence).c_str());
        } // bdtSignalLowerValues
      } // weightPresences
    } // promptnesses
  };

  ProcessMerge(true);
  ProcessMerge(false);

  fileOut->Close();

  return 0;
}
