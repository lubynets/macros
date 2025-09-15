//
// Created by oleksii on 12.09.2025.
//

#include "HelperGeneral.hpp"
#include "HelperMath.hpp"

#include <TAxis.h>
#include <THnSparse.h>
#include <TH1.h>

#include <iostream>

using namespace HelperGeneral;
using namespace HelperMath;

const std::string lifetimeAxisTitle = "T_{proper} (ps)";

const std::vector<float> pTRanges = {0, 2, 5, 8, 12, 20};
const std::string pTAxisTitle = "#it{p}_{T}(#Lambda_{c}^{+}) (GeV/#it{c})";

const std::vector<float> bdtBgUpperValuesVsPt = {0.02, 0.02, 0.02, 0.05, 0.08};
const std::string bgAxisTitle = "BDT bkg score (Lc)";
const std::string npAxisTitle = "BDT non-prompt score (Lc)";
const std::string signalTypeAxisTitle = "candidates type";

std::vector<std::pair<std::string, int>> promptnesses {
  {"prompt", 1},
  {"nonprompt", 2}
};
const std::vector<std::string> weightsPresences{"", "_W"};

std::string GetPtCutName(size_t iPt);

void FillYieldRec(const std::string& fileName, const std::string& filePtWeightName) {
  if(bdtBgUpperValuesVsPt.size() != pTRanges.size() - 1) throw std::runtime_error("bdtUpperValuesVsPt.size() != pTRanges.size() - 1");

  TFile* fileIn = OpenFileWithNullptrCheck(fileName);
  TFile* fileOut = TFile::Open("yield_lifetime_qa_thn.root", "recreate");
  const bool isDoWeight = !filePtWeightName.empty();
  TFile* fileWeight = isDoWeight ? OpenFileWithNullptrCheck(filePtWeightName) : nullptr;
  TH1* histoWeight = isDoWeight ? GetObjectWithNullptrCheck<TH1>(fileWeight, "histoWeight_pT_0_20") : nullptr;

  THnSparse* histoRec = GetObjectWithNullptrCheck<THnSparse>(fileIn, "hf-task-lc/hnLcVarsWithBdt");
  const std::map<std::string, int> axesIndices = MapAxesIndices(histoRec);
  THnSparse* histoRecWeighted = dynamic_cast<THnSparse*>(histoRec->Clone());
  if(isDoWeight) {
    ScaleTHnSparseWithWeight(histoRecWeighted, axesIndices.at(pTAxisTitle), histoWeight);
  }
  
  auto ProcessTHnSparse = [&](THnSparse* histoIn, const std::string& histoNameSuffix="") {
    CheckTAxisForRanges(*histoIn->GetAxis(axesIndices.at(pTAxisTitle)), pTRanges);
    CheckTAxisForRanges(*histoIn->GetAxis(axesIndices.at(bgAxisTitle)), bdtBgUpperValuesVsPt);
    CheckTAxisForRanges(*histoIn->GetAxis(axesIndices.at(signalTypeAxisTitle)), {1, 2, 3});

    std::vector<float> bdtSignalLowerValues;
    for (int iB = 0; iB <= 99; iB++) {
      bdtSignalLowerValues.emplace_back(0.01 * iB);
    }

    for(size_t iPt=0, nPts=pTRanges.size()-1; iPt<nPts; ++iPt) {
      SetTHnSparseAxisRanges(histoIn, axesIndices.at(pTAxisTitle), pTRanges.at(iPt), pTRanges.at(iPt + 1));
      SetTHnSparseAxisRanges(histoIn, axesIndices.at(bgAxisTitle), 0., bdtBgUpperValuesVsPt.at(iPt));
      for(const auto& promptness : promptnesses) {
        const std::string dirName = "rec/" + promptness.first + "/" + GetPtCutName(iPt);
        SetTHnSparseAxisRanges(histoIn, axesIndices.at(signalTypeAxisTitle), static_cast<float>(promptness.second), static_cast<float>(promptness.second)+1.f);
        for(const auto& bsc : bdtSignalLowerValues) {
          const std::string histoName = "hT_NPgt" + to_string_with_precision(bsc, 2) + histoNameSuffix;
          SetTHnSparseAxisRanges(histoIn, axesIndices.at(npAxisTitle), bsc, 1.);
          TH1* histoYield = histoIn->Projection(axesIndices.at(lifetimeAxisTitle));
          histoYield->SetDirectory(nullptr);
          CD(fileOut, dirName);
          histoYield->Write(histoName.c_str());
          SetTHnSparseAxisRanges(histoIn, axesIndices.at(npAxisTitle));
        } // bdtBgUpperValuesVsPt
        SetTHnSparseAxisRanges(histoIn, axesIndices.at(signalTypeAxisTitle));
      } // promptnesses
      SetTHnSparseAxisRanges(histoIn, axesIndices.at(pTAxisTitle));
      SetTHnSparseAxisRanges(histoIn, axesIndices.at(bgAxisTitle));
    } // pTRanges
  };
  
  ProcessTHnSparse(histoRec);
  if(isDoWeight) ProcessTHnSparse(histoRecWeighted, "_W");

  fileOut->Close();
  fileIn->Close();
  if(isDoWeight) fileWeight->Close();
}

void FillYieldGen(const std::string& fileName, const std::string& filePtWeightName) {
  TFile* fileIn = OpenFileWithNullptrCheck(fileName);
  TFile* fileOut = TFile::Open("yield_lifetime_qa_thn.root", "update");
  const bool isDoWeight = !filePtWeightName.empty();
  TFile* fileWeight = isDoWeight ? OpenFileWithNullptrCheck(filePtWeightName) : nullptr;
  TH1* histoWeight = isDoWeight ? GetObjectWithNullptrCheck<TH1>(fileWeight, "histoWeight_pT_0_20") : nullptr;

  THnSparse* histoSim = GetObjectWithNullptrCheck<THnSparse>(fileIn, "hf-task-lc/hnLcVarsGen");
  const std::map<std::string, int> axesIndices = MapAxesIndices(histoSim);
  THnSparse* histoSimWeighted = dynamic_cast<THnSparse*>(histoSim->Clone());
  if(isDoWeight) {
    ScaleTHnSparseWithWeight(histoSimWeighted, axesIndices.at(pTAxisTitle), histoWeight);
  }

  auto ProcessTHnSparse = [&](THnSparse* histoIn, const std::string& histoNameSuffix="") {
    CheckTAxisForRanges(*histoIn->GetAxis(axesIndices.at(pTAxisTitle)), pTRanges);
    CheckTAxisForRanges(*histoIn->GetAxis(axesIndices.at(signalTypeAxisTitle)), {1, 2, 3});

    for(size_t iPt=0, nPts=pTRanges.size()-1; iPt<nPts; ++iPt) {
      SetTHnSparseAxisRanges(histoIn, axesIndices.at(pTAxisTitle), pTRanges.at(iPt), pTRanges.at(iPt + 1));
      for(const auto& promptness : promptnesses) {
        const std::string dirName = "gen/" + promptness.first + "/" + GetPtCutName(iPt);
        SetTHnSparseAxisRanges(histoIn, axesIndices.at(signalTypeAxisTitle), static_cast<float>(promptness.second), static_cast<float>(promptness.second)+1.f);
        const std::string histoName = "hT" + histoNameSuffix;
        TH1* histoYield = histoIn->Projection(axesIndices.at(lifetimeAxisTitle));
        histoYield->SetDirectory(nullptr);
        CD(fileOut, dirName);
        histoYield->Write(histoName.c_str());
        SetTHnSparseAxisRanges(histoIn, axesIndices.at(signalTypeAxisTitle));
      } // promptnesses
      SetTHnSparseAxisRanges(histoIn, axesIndices.at(pTAxisTitle));
    } // pTRanges
  };

  ProcessTHnSparse(histoSim);
  if(isDoWeight) ProcessTHnSparse(histoSimWeighted, "_W");

  fileOut->Close();
  fileIn->Close();
  if(isDoWeight) fileWeight->Close();
}

std::string GetPtCutName(size_t iPt) {
  std::pair<size_t, size_t> iPTMinMax = (iPt == pTRanges.size()-1) ? std::pair<size_t, size_t>{0, pTRanges.size()-1} : std::pair<size_t, size_t>{iPt, iPt+1};
  return "pT_" + to_string_with_precision(pTRanges.at(iPTMinMax.first), 0) + "_" + to_string_with_precision(pTRanges.at(iPTMinMax.second), 0);
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./yield_lifetime_qa_thn fileName (filePtWeightName)" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileName = argv[1];
  const std::string filePtWeightName = argc > 2 ? argv[2] : "";

  FillYieldRec(fileName, filePtWeightName);
  FillYieldGen(fileName, filePtWeightName);

  std::vector<std::string> pTCutNames;
  for(size_t iPt=0, nPts=pTRanges.size()-1; iPt<nPts; ++iPt) {
    pTCutNames.emplace_back(GetPtCutName(iPt));
  }
  TFile* fileOut = OpenFileWithNullptrCheck("yield_lifetime_qa_thn.root", "update");
  for(const auto& promptness : promptnesses) {
    for(const auto& weightPresence : weightsPresences) {
      if(promptness.first == "nonprompt" && weightPresence == "_W") continue;
      std::vector<std::string> histoGenNames;
      histoGenNames.reserve(pTCutNames.size());
      for (const auto& ptcn : pTCutNames) {
        histoGenNames.emplace_back("gen/" + promptness.first + "/" + ptcn + "/hT" + weightPresence);
      }
      TH1* histoGenMerged = MergeHistograms(fileOut, histoGenNames);
      CD(fileOut, "gen/" + promptness.first + "/" + GetPtCutName(pTRanges.size() - 1));
      histoGenMerged->Write(("hT" + weightPresence).c_str());
    } // weightPresences

    std::vector<float> bdtSignalLowerValues;
    for (int iB = 0; iB <= 99; iB++) {
      bdtSignalLowerValues.emplace_back(0.01 * iB);
    }
    for (const auto& bslv : bdtSignalLowerValues) {
      for(const auto& weightPresence : weightsPresences) {
        if(promptness.first == "nonprompt" && weightPresence == "_W") continue;
        std::vector<std::string> histoRecNames;
        histoRecNames.reserve(pTCutNames.size());
        for (const auto& ptcn : pTCutNames) {
          histoRecNames.emplace_back("rec/" + promptness.first + "/" + ptcn + "/hT_NPgt" + to_string_with_precision(bslv, 2) + weightPresence);
        }
        TH1* histoRecMerged = MergeHistograms(fileOut, histoRecNames);
        CD(fileOut, "rec/" + promptness.first + "/" + GetPtCutName(pTRanges.size() - 1));
        histoRecMerged->Write(("hT_NPgt" + to_string_with_precision(bslv, 2) + weightPresence).c_str());
      } // weightPresences
    } // bdtSignalLowerValues
  } // promptnesses
  fileOut->Close();

  return 0;
}