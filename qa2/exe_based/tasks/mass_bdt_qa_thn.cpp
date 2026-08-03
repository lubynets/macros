//
// Created by oleksii on 26.08.2025.
//

#include "HelperGeneral.hpp"
#include "HelperMath.hpp"

#include <TAxis.h>
#include <THnSparse.h>
#include <TH1.h>

#include <iostream>
#include <string>
#include <string_view>
#include <utility>

using namespace HelperGeneral;

const std::vector<double> lifetimeRanges = {0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8};
const std::string lifetimeAxisTitle = "T_{proper} (ps)";

std::vector<double> pTRanges = {3, 4, 5, 8, 12, 20};
const std::string pTAxisTitle = "#it{p}_{T}(#Lambda_{c}^{+}) (GeV/#it{c})";

const std::vector<double> bdtBgUpperValuesVsPt = {0.02, 0.02, 0.02, 0.04, 0.08}; // standard
// const std::vector<double> bdtBgUpperValuesVsPt = {0.01, 0.01, 0.01, 0.02, 0.04}; // tight
// const std::vector<double> bdtBgUpperValuesVsPt = {0.03, 0.03, 0.03, 0.06, 0.12}; // loose
// const std::vector<double> bdtBgUpperValuesVsPt = {0.015, 0.015, 0.015, 0.03, 0.06}; // semitight
// const std::vector<double> bdtBgUpperValuesVsPt = {0.05, 0.05, 0.05, 0.08, 0.16}; // veryloose
// const std::vector<double> bdtBgUpperValuesVsPt = {0.025, 0.025, 0.025, 0.05, 0.10}; // semiloose

const std::string bgAxisTitle = "BDT bkg score (Lc)";
const std::string npAxisTitle = "BDT non-prompt score (Lc)";
const std::string massAxisTitle = "inv. mass (p K #pi) (GeV/#it{c}^{2})";

const std::string bdtScanAxisTitle = npAxisTitle;
const std::string bdtScanDir = "gt";

// const std::string bdtScanAxisTitle = bgAxisTitle;
// const std::string bdtScanDir = "lt";

constexpr int Verobsity{0};

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

  if(bdtScanAxisTitle == bgAxisTitle && modeRun != RunOnly) throw std::runtime_error("bdtScanAxisTitle == bgAxisTitle && modeRun != RunOnly");

  TFile* fileIn = modeRun != MergeOnly ? OpenFileWithNullptrCheck(fileName) : nullptr;
  const std::string fileOutName = modeRun != MergeOnly ? "mass_bdt_qa_thn.root" : fileNameIn;
  const std::string fileOutOption = modeRun != MergeOnly ? "recreate" : "update";
  TFile* fileOut = TFile::Open(fileOutName.c_str(), fileOutOption.c_str());

  std::string bdtScanShortCut;
  if     (bdtScanAxisTitle == npAxisTitle) bdtScanShortCut = "NP";
  else if(bdtScanAxisTitle == bgAxisTitle) bdtScanShortCut = "BG";

  THnSparse* histoIn = modeRun != MergeOnly ? GetObjectWithNullptrCheck<THnSparse>(fileIn, "hf-task-lc/hnLcVarsWithBdt") : nullptr;

  const std::map<std::string_view, int> axesIndices = modeRun != MergeOnly ? MapTHnSparseAxesIndices(histoIn) : std::map<std::string_view, int>{};

  if(bdtBgUpperValuesVsPt.size() != pTRanges.size() - 1) throw std::runtime_error("bdtUpperValuesVsPt.size() != pTRanges.size() - 1");
  if(bdtScanDir != "gt" && bdtScanDir != "lt") throw std::runtime_error("bdtScanDir != \"gt\" && bdtScanDir != \"lt\"");

  std::vector<double> bdtScanValues{0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90};
//   std::vector<double> bdtScanValues{};
//   for (int iB = 0; iB <= 99; iB++) {
//     bdtScanValues.emplace_back(0.01 * iB);
//   }

  std::vector<double> bdtNpUpperValues{1.00, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60};

  if(modeRun != MergeOnly) {
    CheckTAxisForRanges(*histoIn->GetAxis(axesIndices.at(pTAxisTitle)), pTRanges);
    CheckTAxisForRanges(*histoIn->GetAxis(axesIndices.at(bgAxisTitle)), bdtBgUpperValuesVsPt);
    CheckTAxisForRanges(*histoIn->GetAxis(axesIndices.at(lifetimeAxisTitle)), lifetimeRanges);
    CheckTAxisForRanges(*histoIn->GetAxis(axesIndices.at(bdtScanAxisTitle)), bdtScanValues);
    CheckTAxisForRanges(*histoIn->GetAxis(axesIndices.at(bdtScanAxisTitle)), bdtNpUpperValues);
  }

  std::vector<std::string> pTCutNames, tCutNames;
  for(size_t iPt=0, nPts=pTRanges.size()-1; iPt<nPts; ++iPt) {
    pTCutNames.emplace_back(GetPtCutName(iPt));
  }
  for(size_t iT=0, nTs=lifetimeRanges.size()-1; iT<nTs; ++iT) {
    tCutNames.emplace_back("T_" + to_string_with_precision(lifetimeRanges.at(iT), 2) + "_" + to_string_with_precision(lifetimeRanges.at(iT+1), 2));
  }

  for(size_t iPt=0, nPts=pTRanges.size()-1; iPt<nPts && modeRun != MergeOnly; ++iPt) {
    if(Verobsity >=1) std::cout << "\nProcessing iPt = " << iPt << "\n";
    SetTHnSparseAxisRanges(histoIn, axesIndices.at(pTAxisTitle), pTRanges.at(iPt), pTRanges.at(iPt + 1));
    SetTHnSparseAxisRanges(histoIn, axesIndices.at(bgAxisTitle), 0., bdtBgUpperValuesVsPt.at(iPt));
    for(size_t iT=0, nTs=lifetimeRanges.size()-1; iT<nTs; ++iT) {
      if(Verobsity >= 2) std::cout << "Processing iT = " << iT << "\n";
      SetTHnSparseAxisRanges(histoIn, axesIndices.at(lifetimeAxisTitle), lifetimeRanges.at(iT), lifetimeRanges.at(iT + 1));
      if(Verobsity >= 3) std::cout << "Processing bdtNpUpper = ";
      const std::string dirName = pTCutNames.at(iPt) + "/" + tCutNames.at(iT);
      for(const auto& bdtNpUpper : bdtNpUpperValues) {
        if(Verobsity >= 3) std::cout << bdtNpUpper << " ";
        if(Verobsity >= 3) std::cout << "Processing bdtScan = ";
        for (const auto& bdtScan: bdtScanValues) {
          if(Verobsity >= 3) std::cout << bdtScan << " ";
          if(bdtScanAxisTitle == bgAxisTitle && bdtScan > bdtBgUpperValuesVsPt.at(iPt)+0.001) continue;
          if(bdtScanAxisTitle == npAxisTitle && bdtNpUpper <= bdtScan) continue;
          auto [bdtFrom, bdtTo] = bdtScanDir == "gt" ? std::make_pair(bdtScan, 1.) : std::make_pair(0., bdtScan);
          if(bdtScanAxisTitle == npAxisTitle) bdtTo = bdtNpUpper;
          SetTHnSparseAxisRanges(histoIn, axesIndices.at(bdtScanAxisTitle), bdtFrom, bdtTo);
          TH1D* histoMass = histoIn->Projection(axesIndices.at(massAxisTitle));
          histoMass->SetDirectory(nullptr);
          std::string histoName = "hM_" + bdtScanShortCut + bdtScanDir + to_string_with_precision(bdtScan, 2) + "_NPlt" + to_string_with_precision(bdtNpUpper, 2);
          CD(fileOut, dirName);
          histoMass->Write(histoName.c_str());
          SetTHnSparseAxisRanges(histoIn, axesIndices.at(bdtScanAxisTitle));
        } // bdtScanValues
        if(Verobsity >= 3) std::cout << "\n";
      } // bdtNpUpperValues
      if(Verobsity >= 3) std::cout << "\n";
      SetTHnSparseAxisRanges(histoIn, axesIndices.at(lifetimeAxisTitle));
      if(Verobsity >= 2) std::cout << "\n";
    } // lifetimeRanges
    SetTHnSparseAxisRanges(histoIn, axesIndices.at(pTAxisTitle));
    SetTHnSparseAxisRanges(histoIn, axesIndices.at(bgAxisTitle));
  } // pTRanges

  if(modeRun != RunOnly) {
    const double mergePtFrom{3.};
    const double mergePtTo{20};

    const int skipFirstNBins = std::distance(pTRanges.begin(), std::find(pTRanges.begin(), pTRanges.end(), mergePtFrom));
    const int skipLastNBins = std::distance(std::find(pTRanges.begin(), pTRanges.end(), mergePtTo), pTRanges.end()) - 1;

    if(skipFirstNBins == pTRanges.size() || skipLastNBins == -1) throw std::runtime_error("MassBdtQaThn(): mergePtFrom or mergePtTo does not correspond to any of pTRanges");

    pTRanges.erase(pTRanges.end()-skipLastNBins, pTRanges.end());
    pTCutNames.erase(pTCutNames.end()-skipLastNBins, pTCutNames.end());
    pTRanges.erase(pTRanges.begin(), pTRanges.begin()+skipFirstNBins);
    pTCutNames.erase(pTCutNames.begin(), pTCutNames.begin()+skipFirstNBins);

    for (const auto& tcn : tCutNames) {
      for(const auto& bdtNpUpper : bdtNpUpperValues) {
        for (const auto& bslv : bdtScanValues) {
          if(bdtNpUpper <= bslv) continue;
          std::vector<std::string> histoNames;
          histoNames.reserve(pTCutNames.size());
          for (const auto& ptcn : pTCutNames) {
            histoNames.emplace_back(ptcn + "/" + tcn + "/hM_" + bdtScanShortCut + bdtScanDir + to_string_with_precision(bslv, 2) + "_NPlt" + to_string_with_precision(bdtNpUpper, 2));
          }
          TH1* histoMerged = HelperMath::MergeHistograms(fileOut, histoNames);
          HelperGeneral::CD(fileOut, GetPtCutName(pTRanges.size()-1) + "/" + tcn);
          histoMerged->Write(("hM_" + bdtScanShortCut + bdtScanDir + to_string_with_precision(bslv, 2) + "_NPlt" + to_string_with_precision(bdtNpUpper, 2)).c_str());
        } // bdtScanValues
      } // bdtNpUpperValues
    } // TCuts
  } // modeRun != RunOnly

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
  const int modeRun = argc > 2 ? std::stoi(argv[2]) : RunOnly;
  if(modeRun < 0 || modeRun >= NModeRuns) throw std::runtime_error("modeRun < 0 || modeRun >= NModeRuns");

  MassBdtQaThn(fileNameIn, modeRun);

  return 0;
}
