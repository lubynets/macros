//
// Created by oleksii on 26.08.2025.
//

#include "HelperGeneral.hpp"

#include <TAxis.h>
#include <THnSparse.h>

#include <iostream>

using namespace HelperGeneral;

const std::vector<float> pTRanges = {0, 2, 5, 8, 12, 20};
const std::string pTAxisTitle = "#it{p}_{T}(#Lambda_{c}^{+}) (GeV/#it{c})";

const std::vector<float> bdtBgUpperValuesVsPt = {0.02, 0.02, 0.02, 0.05, 0.08};
const std::string bgAxisTitle = "BDT bkg score (Lc)";
const std::string npAxisTitle = "BDT non-prompt score (Lc)";

std::map<std::string, int> MapAxesIndices(THnSparse* histo);
void CheckAxisForRanges(const TAxis& axis, const std::vector<float>& ranges);

void MassBdtQaThn(const std::string& fileName) {
  LoadMacro("styles/mc_qa2.style.cc");

  TFile* fileIn = OpenFileWithNullptrCheck(fileName);

  THnSparse* histoIn = GetObjectWithNullptrCheck<THnSparse>(fileIn, "hf-task-lc/hnLcVarsWithBdt");
  const std::map<std::string, int> axesIndices = MapAxesIndices(histoIn);

  CheckAxisForRanges(*histoIn->GetAxis(axesIndices.at(pTAxisTitle)), pTRanges);
  CheckAxisForRanges(*histoIn->GetAxis(axesIndices.at(bgAxisTitle)), bdtBgUpperValuesVsPt);
  if(bdtBgUpperValuesVsPt.size() != pTRanges.size() - 1) throw std::runtime_error("bdtUpperValuesVsPt.size() != pTRanges.size() - 1");

  fileIn->Close();
}

std::map<std::string, int> MapAxesIndices(THnSparse* histo) {
  std::map<std::string, int> result;
  const short nDims = histo->GetNdimensions();
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