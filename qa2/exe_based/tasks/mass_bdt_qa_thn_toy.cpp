//
// Created by oleksii on 29.08.2025.
//

#include "HelperGeneral.hpp"
#include "HelperMath.hpp"

#include <TAxis.h>
#include <THnSparse.h>
#include <TH1.h>

#include <iostream>

using namespace HelperGeneral;

const std::string lifetimeAxisTitle = "T_{proper} (ps)";
const std::string pTAxisTitle = "#it{p}_{T}(#Lambda_{c}^{+}) (GeV/#it{c})";
const std::string bgAxisTitle = "BDT bkg score (Lc)";
const std::string npAxisTitle = "BDT non-prompt score (Lc)";
const std::string massAxisTitle = "inv. mass (p K #pi) (GeV/#it{c}^{2})";

std::map<std::string, int> MapAxesIndices(const THnSparse* histo);
void CheckAxisForRanges(const TAxis& axis, const std::vector<float>& ranges);
void SetRanges(THnSparse* histo, int axisNum, float lo=-999., float hi=-999.);

void MassBdtQaThnToy(const std::string& fileName) {
  TFile* fileIn = OpenFileWithNullptrCheck(fileName);
  TFile* fileOut = TFile::Open("mass_bdt_qa_thn_toy.root", "recreate");

  THnSparse* histoIn = GetObjectWithNullptrCheck<THnSparse>(fileIn, "hf-task-lc/hnLcVarsWithBdt");
  const std::map<std::string, int> axesIndices = MapAxesIndices(histoIn);

  TH1* hPt = histoIn->Projection(axesIndices.at(pTAxisTitle));
  TH1* hMass = histoIn->Projection(axesIndices.at(massAxisTitle));
  TH1* hBg = histoIn->Projection(axesIndices.at(bgAxisTitle));
  TH1* hNp = histoIn->Projection(axesIndices.at(npAxisTitle));
  fileOut->cd();
  hPt->Write("hPt");
  hMass->Write("hMass");
  hBg->Write("hBg");
  hNp->Write("hNp");

  for(int iPt=1; iPt<10; ++iPt) {
    const float lo = 0.1 * iPt;
    const float hi = 0.1 * (iPt + 1);

    SetRanges(histoIn, axesIndices.at(pTAxisTitle), lo, hi);
    TH1* histoMass = histoIn->Projection(axesIndices.at(massAxisTitle));
    histoMass->SetDirectory(nullptr);
    fileOut->cd();
    const std::string histoName = "hMass_" + to_string_with_precision(lo, 1) + "_" + to_string_with_precision(hi, 1);
    histoMass->Write(histoName.c_str());
  }

  fileOut->Close();
  fileIn->Close();
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
  std::cout << "binLo = " << binLo << "; binHi = " << binHi << "\n";
}

std::map<std::string, int> MapAxesIndices(const THnSparse* histo) {
  std::map<std::string, int> result;
  const int nDims = histo->GetNdimensions();
  for(int iDim=0; iDim<nDims; ++iDim) {
    result.insert({histo->GetAxis(iDim)->GetTitle(), iDim});
  }
  return result;
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./mass_bdt_qa_thn_toy fileName" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileName = argv[1];

  MassBdtQaThnToy(fileName);

  return 0;
}
