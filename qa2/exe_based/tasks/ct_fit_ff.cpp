//
// Created by oleksii on 14.09.2026
//
#include "ForwardFoldedExpo.hpp"
#include "HelperCtFit.hpp"
#include "HelperGeneral.hpp"

#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

constexpr HelperCtFit::RunModes RunMode{
//   HelperCtFit::RunModes::MeanFitOnly
//   HelperCtFit::RunModes::AllWoOne
  HelperCtFit::RunModes::AllPossible
};

void ExtendHistoWithEmptyBinsLeft(TH1*& histo, const std::vector<double>& additionalBinEdges);

[[deprecated]] void zeroNonDiagonalBins(TH2* histo) {
  const int nBins = histo->GetNbinsX();
  for(int iBin=1; iBin<=nBins; ++iBin) {
    for(int jBin=1; jBin<=nBins; ++jBin) {
      if(iBin != jBin) histo->SetBinContent(iBin, jBin, 0.);
    }
  }
}

void ct_fit_ff(const std::string& fileNameYield, const std::string& fileNameResponseMatrix, const std::string& histoNameResponseMatrix) {
  TFile* fileYield = HelperGeneral::OpenFileWithNullptrCheck(fileNameYield);
  TH1* histoYield = HelperGeneral::GetObjectWithNullptrCheck<TH1>(fileYield, "hCorrYieldsPrompt");
  ExtendHistoWithEmptyBinsLeft(histoYield, {0., 0.2});

  histoYield->SetMarkerColor(kBlue);
  histoYield->SetLineColor(kBlue);
  histoYield->SetMarkerStyle(kFullSquare);
  histoYield->SetMarkerSize(1.6);
  histoYield->GetYaxis()->SetTitle("semi-corrected yield prompt");

  TFile* fileRespMatrix = HelperGeneral::OpenFileWithNullptrCheck(fileNameResponseMatrix);
  TH2* histoRespMatrix = HelperGeneral::GetObjectWithNullptrCheck<TH2>(fileRespMatrix, histoNameResponseMatrix);
//   zeroNonDiagonalBins(histoRespMatrix);

  HelperGeneral::CheckHistogramsForAxisIdentity<TH2, TH2>(histoRespMatrix, nullptr, "XY");
  HelperGeneral::CheckHistogramsForAxisIdentity(histoRespMatrix, histoYield, "X");

  const double lo = 0.2;
  const double hi = 1.8;

  // ------------- expo parameters estimate --------------------------------------------
  TH1* histEff = histoRespMatrix->ProjectionX(); // for expo parameters estimate only
  TH1* histPreliminaryCorrectedYield = dynamic_cast<TH1*>(histoYield->Clone());
  histPreliminaryCorrectedYield->SetDirectory(nullptr);
  histPreliminaryCorrectedYield->Divide(histEff);
  const double preliminaryYield = histPreliminaryCorrectedYield->Integral();
  const double preliminaryA = preliminaryYield / HelperGeneral::LifetimeLambdaC / (std::exp(-lo / HelperGeneral::LifetimeLambdaC) - std::exp(-hi / HelperGeneral::LifetimeLambdaC));
  // -----------------------------------------------------------------------------------

  ForwardFoldedExpo ffe{};
  ffe.SetResponseMatrix(histoRespMatrix);
  ffe.Init();

  TF1* fitFunc = new TF1("fitFunc", ffe, 0., 2., 2);
  fitFunc->SetParameters(preliminaryA, HelperGeneral::LifetimeLambdaC);
  fitFunc->SetNpx(1000);

  histoYield->Fit(fitFunc, "", "", lo, hi);

  histoYield->SaveAs("h1.root");

  fileRespMatrix->Close();
  fileYield->Close();
}

void ExtendHistoWithEmptyBinsLeft(TH1*& histo, const std::vector<double>& additionalBinEdges) {
  if(additionalBinEdges.size() < 2) {
    throw std::runtime_error("ExtendHistoWithEmptyBinsLeft() - need at least 2 bin edges");
  }
  if(!HelperGeneral::EqualFloating(additionalBinEdges.back(), histo->GetXaxis()->GetXmin())) {
    throw std::runtime_error("ExtendHistoWithEmptyBinsLeft() - additionalBinEdges.back() != histo->GetXaxis()->GetXmin()");
  }
  if(!std::is_sorted(additionalBinEdges.begin(), additionalBinEdges.end())) {
    throw std::runtime_error("ExtendHistoWithEmptyBinsLeft() - additionalBinEdges is not sorted");
  }
  const int nAddBins = additionalBinEdges.size() - 1;
  std::vector<double> newHistoEdges{additionalBinEdges};
  for(int iBin=1, nBins=histo->GetNbinsX(); iBin<=nBins; ++iBin) {
    newHistoEdges.push_back(histo->GetBinLowEdge(iBin+1));
  }
  TH1* hResult = new TH1D(histo->GetName(), histo->GetTitle(), newHistoEdges.size()-1, newHistoEdges.data());
  hResult->GetXaxis()->SetTitle(histo->GetXaxis()->GetTitle());
  hResult->GetYaxis()->SetTitle(histo->GetYaxis()->GetTitle());
  const bool hasSumw2 = histo->GetSumw2N() > 0;
  if(hasSumw2) hResult->Sumw2();
  for(int iBin=1, nBins=histo->GetNbinsX(); iBin<=nBins; ++iBin) {
    hResult->SetBinContent(iBin+nAddBins, histo->GetBinContent(iBin));
    if(hasSumw2) hResult->SetBinError(iBin+nAddBins, histo->GetBinError(iBin));
  }
  histo = hResult;
}

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./ct_fit_ff fileNameYield fileNameResponseMatrix histoNameResponseMatrix" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileNameYield = argv[1];
  const std::string fileNameResponseMatrix = argv[2];
  const std::string histoNameResponseMatrix = argv[3];

  ct_fit_ff(fileNameYield, fileNameResponseMatrix, histoNameResponseMatrix);

  return 0;
}
