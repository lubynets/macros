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

const std::string lifetimeAxisTitle = "T_{proper} (ps)";

const std::vector<float> pTRanges = {0, 2, 5, 8, 12, 20};
const std::string pTAxisTitle = "#it{p}_{T}(#Lambda_{c}^{+}) (GeV/#it{c})";

const std::vector<float> bdtBgUpperValuesVsPt = {0.02, 0.02, 0.02, 0.05, 0.08};
const std::string bgAxisTitle = "BDT bkg score (Lc)";
const std::string npAxisTitle = "BDT non-prompt score (Lc)";
const std::string massAxisTitle = "inv. mass (p K #pi) (GeV/#it{c}^{2})";

std::string GetPtCutName(size_t iPt);

void FillYieldRec(const std::string& fileName) {
  TFile* fileIn = OpenFileWithNullptrCheck(fileName);
  TFile* fileOut = TFile::Open("yield_lifetime_qa_thn.root", "recreate");

  THnSparse* histoIn = GetObjectWithNullptrCheck<THnSparse>(fileIn, "hf-task-lc/hnLcVarsWithBdt");
  const std::map<std::string, int> axesIndices = MapAxesIndices(histoIn);

  CheckTAxisForRanges(*histoIn->GetAxis(axesIndices.at(pTAxisTitle)), pTRanges);
  CheckTAxisForRanges(*histoIn->GetAxis(axesIndices.at(bgAxisTitle)), bdtBgUpperValuesVsPt);
  if(bdtBgUpperValuesVsPt.size() != pTRanges.size() - 1) throw std::runtime_error("bdtUpperValuesVsPt.size() != pTRanges.size() - 1");

  std::vector<float> bdtSignalLowerValues;
  for (int iB = 0; iB <= 99; iB++) {
    bdtSignalLowerValues.emplace_back(0.01 * iB);
  }
  std::vector<std::string> pTCutNames, tCutNames;
  for(size_t iPt=0, nPts=pTRanges.size()-1; iPt<nPts; ++iPt) {
    pTCutNames.emplace_back(GetPtCutName(iPt));
  }




  fileOut->Close();
  fileIn->Close();
}

void FillYieldGen(const std::string& fileName) {

}

std::string GetPtCutName(size_t iPt) {
  std::pair<size_t, size_t> iPTMinMax = (iPt == pTRanges.size()-1) ? std::pair<size_t, size_t>{0, pTRanges.size()-1} : std::pair<size_t, size_t>{iPt, iPt+1};
  return "pT_" + to_string_with_precision(pTRanges.at(iPTMinMax.first), 0) + "_" + to_string_with_precision(pTRanges.at(iPTMinMax.second), 0);
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./yield_lifetime_qa_thn fileName" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileName = argv[1];

  FillYieldRec(fileName);
  FillYieldGen(fileName);

  return 0;
}