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

const std::vector<float> lifetimeRanges = {0.2, 0.35, 0.5, 0.7, 0.9, 1.6};
const std::string lifetimeAxisTitle = "T_{proper} (ps)";

const std::vector<float> pTRanges = {0, 2, 5, 8, 12, 20};
const std::string pTAxisTitle = "#it{p}_{T}(#Lambda_{c}^{+}) (GeV/#it{c})";

const std::vector<float> bdtBgUpperValuesVsPt = {0.02, 0.02, 0.02, 0.05, 0.08};
const std::string bgAxisTitle = "BDT bkg score (Lc)";
const std::string npAxisTitle = "BDT non-prompt score (Lc)";
const std::string massAxisTitle = "inv. mass (p K #pi) (GeV/#it{c}^{2})";

std::map<std::string, int> MapAxesIndices(const THnSparse* histo);
void CheckAxisForRanges(const TAxis& axis, const std::vector<float>& ranges);
std::string GetPtCutName(size_t iPt);
void SetRanges(THnSparse* histo, int axisNum, float lo=-999., float hi=-999.);

void MassBdtQaThn(const std::string& fileName) {
//  LoadMacro("styles/mc_qa2.style.cc");

  TFile* fileIn = OpenFileWithNullptrCheck(fileName);
  TFile* fileOut = TFile::Open("mass_bdt_qa_thn.root", "recreate");

  THnSparse* histoIn = GetObjectWithNullptrCheck<THnSparse>(fileIn, "hf-task-lc/hnLcVarsWithBdt");
  const std::map<std::string, int> axesIndices = MapAxesIndices(histoIn);

  CheckAxisForRanges(*histoIn->GetAxis(axesIndices.at(pTAxisTitle)), pTRanges);
  CheckAxisForRanges(*histoIn->GetAxis(axesIndices.at(bgAxisTitle)), bdtBgUpperValuesVsPt);
  CheckAxisForRanges(*histoIn->GetAxis(axesIndices.at(lifetimeAxisTitle)), lifetimeRanges);
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

  for(size_t iPt=0, nPts=pTRanges.size()-1; iPt<nPts; ++iPt) {
    SetRanges(histoIn, axesIndices.at(pTAxisTitle), pTRanges.at(iPt), pTRanges.at(iPt + 1));
    SetRanges(histoIn, axesIndices.at(bgAxisTitle), 0., bdtBgUpperValuesVsPt.at(iPt));
    for(size_t iT=0, nTs=lifetimeRanges.size()-1; iT<nTs; ++iT) {
      const std::string dirName = pTCutNames.at(iPt) + "/" + tCutNames.at(iT);
      SetRanges(histoIn, axesIndices.at(lifetimeAxisTitle), lifetimeRanges.at(iT), lifetimeRanges.at(iT+1));
      for (const auto& bdtSig: bdtSignalLowerValues) {
        SetRanges(histoIn, axesIndices.at(npAxisTitle), bdtSig, 1.);
        TH1D* histoMass = histoIn->Projection(axesIndices.at(massAxisTitle));
        histoMass->SetDirectory(nullptr);
        const std::string histoName = "hM_NPgt" + to_string_with_precision(bdtSig, 2);
        CD(fileOut, dirName);
        histoMass->Write(histoName.c_str());
        SetRanges(histoIn, axesIndices.at(npAxisTitle));
      } // bdtSignalLowerValues
      SetRanges(histoIn, axesIndices.at(lifetimeAxisTitle));
    } // lifetimeRanges
    SetRanges(histoIn, axesIndices.at(pTAxisTitle));
    SetRanges(histoIn, axesIndices.at(bgAxisTitle));
  } // pTRanges

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

  fileOut->Close();
  fileIn->Close();
}

std::map<std::string, int> MapAxesIndices(const THnSparse* histo) {
  std::map<std::string, int> result;
  const int nDims = histo->GetNdimensions();
  for(int iDim=0; iDim<nDims; ++iDim) {
    result.insert({histo->GetAxis(iDim)->GetTitle(), iDim});
  }
  return result;
}

void CheckAxisForRanges(const TAxis& axis, const std::vector<float>& ranges) {
  const int nBins = axis.GetNbins();
  for(const auto& range : ranges) {
    bool ok{false};
    for(int iBin=1; iBin<=nBins+1; ++iBin) {
      const float edge = axis.GetBinLowEdge(iBin);
      if(std::fabs(edge - range) < 1e-4) {
        ok = true;
        break;
      }
    }
    if(!ok) {
      throw std::runtime_error("CheckAxisForRanges() - the range " + std::to_string(range) + " is missing");
    }
  }
}

std::string GetPtCutName(size_t iPt) {
  std::pair<size_t, size_t> iPTMinMax = (iPt == pTRanges.size()-1) ? std::pair<size_t, size_t>{0, pTRanges.size()-1} : std::pair<size_t, size_t>{iPt, iPt+1};
  return "pT_" + to_string_with_precision(pTRanges.at(iPTMinMax.first), 0) + "_" + to_string_with_precision(pTRanges.at(iPTMinMax.second), 0);
}

void SetRanges(THnSparse* histo, int axisNum, float lo, float hi) {
  constexpr double tolerance = 1e-6;

  if(std::fabs(lo+999)<tolerance && std::fabs(hi+999)<tolerance) {
    histo->GetAxis(axisNum)->SetRange();
    return;
  }

  if(lo >= hi) throw std::runtime_error("SetRanges(): lo >= hi");

  const TAxis* axis = histo->GetAxis(axisNum);
  int binLo{-999}, binHi{-999};
  for(int iBin=1, nBins=axis->GetNbins(); iBin<=nBins; ++iBin) {
    const float binLowEdge = axis->GetBinLowEdge(iBin);
    const float binUpEdge = axis->GetBinUpEdge(iBin);
    if(std::fabs(binLowEdge - lo)<tolerance) binLo = iBin;
    if(std::fabs(binUpEdge - hi)<tolerance) binHi = iBin;
    if(binLo != -999 && binHi != -999) break;
  }
  if(binLo == -999 || binHi == -999) throw std::runtime_error("SetRanges(): binLo == -999 || binHi == -999");
  histo->GetAxis(axisNum)->SetRange(binLo, binHi);
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./mass_bdt_qa_thn fileName" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileName = argv[1];

  MassBdtQaThn(fileName);

  return 0;
}