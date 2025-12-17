//
// Created by oleksii on 26.08.2025.
//

#include "HelperGeneral.hpp"
#include "HelperMath.hpp"

#include <TAxis.h>
#include <THnSparse.h>
#include <TH1.h>

#include <iostream>

using namespace HelperGeneral;

const std::vector<float> lifetimeRanges = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 5.0};
const std::string lifetimeAxisTitle = "T_{proper} (ps)";

std::vector<float> pTRanges = {0, 1, 2, 3, 4, 5, 8, 12, 20};
const std::string pTAxisTitle = "#it{p}_{T}(#Lambda_{c}^{+}) (GeV/#it{c})";

const std::vector<float> bdtBgUpperValuesVsPt = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.1, 0.2};
const std::string bgAxisTitle = "BDT bkg score (Lc)";
const std::string npAxisTitle = "BDT non-prompt score (Lc)";
const std::string massAxisTitle = "inv. mass (p K #pi) (GeV/#it{c}^{2})";

std::string GetPtCutName(size_t iPt);

enum : int {
  RunOnly = 0,
  RunAndMerge,
  MergeOnly,
  NModeRuns
};

void MassBdtQaThn(const std::string& fileNameIn, int modeRun) {
//  LoadMacro("styles/mc_qa2.style.cc");
  const std::string fileName = ReadNthLine(fileNameIn);

  TFile* fileIn = modeRun != MergeOnly ? OpenFileWithNullptrCheck(fileName) : nullptr;
  const std::string fileOutName = modeRun != MergeOnly ? "mass_bdt_qa_thn.root" : fileNameIn;
  const std::string fileOutOption = modeRun != MergeOnly ? "recreate" : "update";
  TFile* fileOut = TFile::Open(fileOutName.c_str(), fileOutOption.c_str());

  THnSparse* histoIn = modeRun != MergeOnly ? GetObjectWithNullptrCheck<THnSparse>(fileIn, "hf-task-lc/hnLcVarsWithBdt") : nullptr;

  const std::map<std::string, int> axesIndices = modeRun != MergeOnly ? MapAxesIndices(histoIn) : std::map<std::string, int>{};
  if(modeRun != MergeOnly) {
    CheckTAxisForRanges(*histoIn->GetAxis(axesIndices.at(pTAxisTitle)), pTRanges);
    CheckTAxisForRanges(*histoIn->GetAxis(axesIndices.at(bgAxisTitle)), bdtBgUpperValuesVsPt);
    CheckTAxisForRanges(*histoIn->GetAxis(axesIndices.at(lifetimeAxisTitle)), lifetimeRanges);
  }

  if(bdtBgUpperValuesVsPt.size() != pTRanges.size() - 1) throw std::runtime_error("bdtUpperValuesVsPt.size() != pTRanges.size() - 1");

  std::vector<float> bdtSignalLowerValues;
  for (int iB = 0; iB <= 99; iB++) {
    bdtSignalLowerValues.emplace_back(0.01 * iB);
  }
  std::vector<std::string> pTCutNames, tCutNames;
  for(size_t iPt=0, nPts=pTRanges.size()-1; iPt<nPts; ++iPt) {
    pTCutNames.emplace_back(GetPtCutName(iPt));
  }
  for(size_t iT=0, nTs=lifetimeRanges.size()-1; iT<nTs; ++iT) {
    tCutNames.emplace_back("T_" + to_string_with_precision(lifetimeRanges.at(iT), 2) + "_" + to_string_with_precision(lifetimeRanges.at(iT+1), 2));
  }

  for(size_t iPt=0, nPts=pTRanges.size()-1; iPt<nPts && modeRun != MergeOnly; ++iPt) {
    std::cout << "\nProcessing iPt = " << iPt << "\n";
    SetTHnSparseAxisRanges(histoIn, axesIndices.at(pTAxisTitle), pTRanges.at(iPt), pTRanges.at(iPt + 1));
    SetTHnSparseAxisRanges(histoIn, axesIndices.at(bgAxisTitle), 0., bdtBgUpperValuesVsPt.at(iPt));
    for(size_t iT=0, nTs=lifetimeRanges.size()-1; iT<nTs; ++iT) {
      std::cout << "Processing iT = " << iT << "\n";
      const std::string dirName = pTCutNames.at(iPt) + "/" + tCutNames.at(iT);
      SetTHnSparseAxisRanges(histoIn, axesIndices.at(lifetimeAxisTitle), lifetimeRanges.at(iT), lifetimeRanges.at(iT + 1));
      std::cout << "Processing bdtSig = ";
      for (const auto& bdtSig: bdtSignalLowerValues) {
        std::cout << bdtSig << " ";
        SetTHnSparseAxisRanges(histoIn, axesIndices.at(npAxisTitle), bdtSig, 1.);
        TH1D* histoMass = histoIn->Projection(axesIndices.at(massAxisTitle));
        histoMass->SetDirectory(nullptr);
        const std::string histoName = "hM_NPgt" + to_string_with_precision(bdtSig, 2);
        CD(fileOut, dirName);
        histoMass->Write(histoName.c_str());
        SetTHnSparseAxisRanges(histoIn, axesIndices.at(npAxisTitle));
      } // bdtSignalLowerValues
      std::cout << "\n";
      SetTHnSparseAxisRanges(histoIn, axesIndices.at(lifetimeAxisTitle));
    } // lifetimeRanges
    SetTHnSparseAxisRanges(histoIn, axesIndices.at(pTAxisTitle));
    SetTHnSparseAxisRanges(histoIn, axesIndices.at(bgAxisTitle));
  } // pTRanges

  if(modeRun != RunOnly) {
    for (const auto& tcn : tCutNames) {
      for (const auto& bslv : bdtSignalLowerValues) {
        std::vector<std::string> histoNames;
        histoNames.reserve(pTCutNames.size());
        for (const auto& ptcn : pTCutNames) {
          histoNames.emplace_back(ptcn + "/" + tcn + "/hM_NPgt" + to_string_with_precision(bslv, 2));
        }
        TH1* histoMerged = HelperMath::MergeHistograms(fileOut, histoNames);
        HelperGeneral::CD(fileOut, GetPtCutName(pTRanges.size()-1) + "/" + tcn);
        histoMerged->Write(("hM_NPgt" + to_string_with_precision(bslv, 2)).c_str());
      } // bdtSignalLowerValues
    } // TCuts
  }

  fileOut->Close();
  if(modeRun != MergeOnly) fileIn->Close();
}

std::string GetPtCutName(size_t iPt) {
  std::pair<size_t, size_t> iPTMinMax = (iPt == pTRanges.size()-1) ? std::pair<size_t, size_t>{0, pTRanges.size()-1} : std::pair<size_t, size_t>{iPt, iPt+1};
  return "pT_" + to_string_with_precision(pTRanges.at(iPTMinMax.first), 0) + "_" + to_string_with_precision(pTRanges.at(iPTMinMax.second), 0);
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./mass_bdt_qa_thn fileNameIn (modeRun=RunOnly=0 [RunAndMerge=1, MergeOnly=2])" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileNameIn = argv[1];
  const int modeRun = argc > 2 ? std::stoi(argv[2]) : RunAndMerge;
  if(modeRun < 0 || modeRun >= NModeRuns) throw std::runtime_error("modeRun < 0 || modeRun >= NModeRuns");

  MassBdtQaThn(fileNameIn, modeRun);

  return 0;
}
