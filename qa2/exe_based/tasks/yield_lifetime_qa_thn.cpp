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

bool gIsDoWeight{false};
std::vector<float> gBdtSignalLowerValues{};

const std::vector<float> pTRanges = {1, 3, 5, 8, 12, 20};
const std::vector<float> bdtBgUpperValuesVsPt = {0.02, 0.02, 0.02, 0.05, 0.08};

const std::string pTAxisTitle = "#it{p}_{T}(#Lambda_{c}^{+}) (GeV/#it{c})";
const std::string bgAxisTitle = "BDT bkg score (Lc)";
const std::string npAxisTitle = "BDT non-prompt score (Lc)";
const std::string signalTypeAxisTitle = "candidates type";

const std::string fileOutName{"yield_lifetime_qa_thn.root"};

std::vector<std::pair<std::string, float>> promptnesses {
  {"prompt", 1},
  {"nonprompt", 2}
};
const std::vector<std::string> weightsPresences{"", "_W"};

std::string GetPtCutName(size_t iPt);

void FillYield(const std::string& fileName, const std::string& filePtWeightName, const bool isRec) {
  TFile* fileIn = OpenFileWithNullptrCheck(fileName);
  const std::string fileOpenOption = isRec ? "recreate" : "update";
  TFile* fileOut = TFile::Open(fileOutName.c_str(), fileOpenOption.c_str());
  TFile* fileWeight = gIsDoWeight ? OpenFileWithNullptrCheck(filePtWeightName) : nullptr;
  TH1* histoWeight = gIsDoWeight ? GetObjectWithNullptrCheck<TH1>(fileWeight, "histoWeight_pT_0_20") : nullptr;

  THnSparse* histoRec = GetObjectWithNullptrCheck<THnSparse>(fileIn, "hf-task-lc/hnLcVarsWithBdt");
  const std::map<std::string, int> axesIndices = MapAxesIndices(histoRec);
  THnSparse* histoRecWeighted = dynamic_cast<THnSparse*>(histoRec->Clone());
  if (gIsDoWeight) {
    ScaleTHnSparseWithWeight(histoRecWeighted, axesIndices.at(pTAxisTitle), histoWeight);
  }
  
  auto ProcessTHnSparse = [&](THnSparse* histoIn, const std::string& histoNameSuffix="", const std::vector<std::pair<std::string, float>>& promptnessesToProcess=promptnesses) {
    CheckTAxisForRanges(*histoIn->GetAxis(axesIndices.at(pTAxisTitle)), pTRanges);
    CheckTAxisForRanges(*histoIn->GetAxis(axesIndices.at(signalTypeAxisTitle)), {1, 2, 3});
    if(isRec) CheckTAxisForRanges(*histoIn->GetAxis(axesIndices.at(bgAxisTitle)), bdtBgUpperValuesVsPt);

    for(size_t iPt=0, nPts=pTRanges.size()-1; iPt<nPts; ++iPt) {
      SetTHnSparseAxisRanges(histoIn, axesIndices.at(pTAxisTitle), pTRanges.at(iPt), pTRanges.at(iPt + 1));
      if(isRec) SetTHnSparseAxisRanges(histoIn, axesIndices.at(bgAxisTitle), 0., bdtBgUpperValuesVsPt.at(iPt));
      for(const auto& promptness : promptnessesToProcess) {
        const std::string dirName = (isRec ? "rec/" : "gen/") + promptness.first + "/" + GetPtCutName(iPt);
        SetTHnSparseAxisRanges(histoIn, axesIndices.at(signalTypeAxisTitle), promptness.second, promptness.second+1.f);
        // for rec - real gBdtSignalLowerValues; for gen - fake 1-element vector for universality reasons
        const auto& bdtSignalLowerValues = isRec ? gBdtSignalLowerValues : std::vector<float>{-999.f};
        for(const auto& bsc : bdtSignalLowerValues) {
          const std::string histoName = isRec ?
                                        "hT_NPgt" + to_string_with_precision(bsc, 2) + histoNameSuffix :
                                        "hT" + histoNameSuffix;
          if(isRec) SetTHnSparseAxisRanges(histoIn, axesIndices.at(npAxisTitle), bsc, 1.);
          TH1* histoYield = histoIn->Projection(axesIndices.at(lifetimeAxisTitle));
          histoYield->SetDirectory(nullptr);
          CD(fileOut, dirName);
          histoYield->Write(histoName.c_str());
          SetTHnSparseAxisRanges(histoIn, axesIndices.at(npAxisTitle));
        } // bdtBgUpperValuesVsPt
        SetTHnSparseAxisRanges(histoIn, axesIndices.at(signalTypeAxisTitle));
      } // promptnessesToProcess
      SetTHnSparseAxisRanges(histoIn, axesIndices.at(pTAxisTitle));
      if(isRec) SetTHnSparseAxisRanges(histoIn, axesIndices.at(bgAxisTitle));
    } // pTRanges
  };
  
  ProcessTHnSparse(histoRec);
  if (gIsDoWeight) ProcessTHnSparse(histoRecWeighted, "_W", {promptnesses.at(0)});

  fileOut->Close();
  fileIn->Close();
  if (gIsDoWeight) fileWeight->Close();
}

std::string GetPtCutName(size_t iPt) {
  std::pair<size_t, size_t> iPTMinMax = (iPt == pTRanges.size()-1) ? std::pair<size_t, size_t>{0, pTRanges.size()-1} : std::pair<size_t, size_t>{iPt, iPt+1};
  return "pT_" + to_string_with_precision(pTRanges.at(iPTMinMax.first), 0) + "_" + to_string_with_precision(pTRanges.at(iPTMinMax.second), 0);
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./yield_lifetime_qa_thn fileNameIn (filePtWeightName)" << std::endl;
    exit(EXIT_FAILURE);
  }
  if(bdtBgUpperValuesVsPt.size() != pTRanges.size() - 1) throw std::runtime_error("bdtUpperValuesVsPt.size() != pTRanges.size() - 1");

  const std::string fileNameIn = argv[1];
  const std::string filePtWeightName = argc > 2 ? argv[2] : "";

  gIsDoWeight = !filePtWeightName.empty();

  const std::string fileName = ReadNthLine(fileNameIn);

  for (int iB = 0; iB <= 99; iB++) {
    gBdtSignalLowerValues.emplace_back(0.01 * iB);
  }

  FillYield(fileName, filePtWeightName, true);
  FillYield(fileName, filePtWeightName, false);

  std::vector<std::string> pTCutNames;
  for(size_t iPt=0, nPts=pTRanges.size()-1; iPt<nPts; ++iPt) {
    pTCutNames.emplace_back(GetPtCutName(iPt));
  }
  TFile* fileOut = OpenFileWithNullptrCheck(fileOutName.c_str(), "update");
  for(const auto& promptness : promptnesses) {
    for(const auto& weightPresence : weightsPresences) {
      if((promptness.first == "nonprompt" || !gIsDoWeight) && weightPresence == "_W") continue;
      std::vector<std::string> histoGenNames;
      histoGenNames.reserve(pTCutNames.size());
      for (const auto& ptcn : pTCutNames) {
        histoGenNames.emplace_back("gen/" + promptness.first + "/" + ptcn + "/hT" + weightPresence);
      }
      TH1* histoGenMerged = MergeHistograms(fileOut, histoGenNames);
      CD(fileOut, "gen/" + promptness.first + "/" + GetPtCutName(pTRanges.size() - 1));
      histoGenMerged->Write(("hT" + weightPresence).c_str());
    } // weightPresences


    for (const auto& bslv : gBdtSignalLowerValues) {
      for(const auto& weightPresence : weightsPresences) {
        if((promptness.first == "nonprompt" || !gIsDoWeight) && weightPresence == "_W") continue;
        std::vector<std::string> histoRecNames;
        histoRecNames.reserve(pTCutNames.size());
        for (const auto& ptcn : pTCutNames) {
          histoRecNames.emplace_back("rec/" + promptness.first + "/" + ptcn + "/hT_NPgt" + to_string_with_precision(bslv, 2) + weightPresence);
        }
        TH1* histoRecMerged = MergeHistograms(fileOut, histoRecNames);
        CD(fileOut, "rec/" + promptness.first + "/" + GetPtCutName(pTRanges.size() - 1));
        histoRecMerged->Write(("hT_NPgt" + to_string_with_precision(bslv, 2) + weightPresence).c_str());
      } // weightPresences
    } // gBdtSignalLowerValues
  } // promptnesses
  fileOut->Close();

  return 0;
}
