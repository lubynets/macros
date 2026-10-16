//
// Created by oleksii on 14.09.2026
//
#include "ForwardFoldedExpo.hpp"
#include "HelperCtFit.hpp"
#include "HelperGeneral.hpp"
#include "HelperMath.hpp"
#include "HelperPlot.hpp"

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

constexpr bool IsSaveCanvasAsRoot{true};

constexpr HelperCtFit::RunModes RunMode{
//   HelperCtFit::RunModes::MeanFitOnly
  HelperCtFit::RunModes::AllWoOne
//   HelperCtFit::RunModes::AllPossible
};

const std::vector<double> additionalBinEdges{0., 0.2};
constexpr double lo = 0.2;
constexpr double hi = 1.8;

TH1* ExtendHistoWithEmptyBinsLeft(TH1* histo, const std::vector<double>& additionalBinEdges);

[[deprecated]] void zeroNonDiagonalBins(TH2* histo) {
  const int nBins = histo->GetNbinsX();
  for(int iBin=1; iBin<=nBins; ++iBin) {
    for(int jBin=1; jBin<=nBins; ++jBin) {
      if(iBin != jBin) histo->SetBinContent(iBin, jBin, 0.);
    }
  }
}

void ct_fit_ff(const std::string& fileNameYield, const std::string& fileNameResponseMatrix, const std::string& histoNameResponseMatrix) {
  HelperGeneral::LoadMacro("styles/mc_qa2.style.cc");

  TFile* fileYield = HelperGeneral::OpenFileWithNullptrCheck(fileNameYield);
  TH1* histoYield = HelperGeneral::GetObjectWithNullptrCheck<TH1>(fileYield, "hCorrYieldsPrompt");
  histoYield->SetDirectory(nullptr);
  histoYield->UseCurrentStyle();
  TH1* histoYieldExt{histoYield};
  if(!additionalBinEdges.empty()) histoYieldExt = ExtendHistoWithEmptyBinsLeft(histoYield, additionalBinEdges);

  histoYieldExt->SetMarkerColor(kBlue);
  histoYieldExt->SetLineColor(kBlue);
  histoYieldExt->SetMarkerStyle(kFullSquare);
  histoYieldExt->SetMarkerSize(1.6);
  histoYieldExt->GetYaxis()->SetTitle("Y_{prompt}^{ref. sel.}");

  TFile* fileRespMatrix = HelperGeneral::OpenFileWithNullptrCheck(fileNameResponseMatrix);
  TH2* histoRespMatrix = HelperGeneral::GetObjectWithNullptrCheck<TH2>(fileRespMatrix, histoNameResponseMatrix);
//   zeroNonDiagonalBins(histoRespMatrix);

  HelperGeneral::CheckHistogramsForAxisIdentity<TH2, TH2>(histoRespMatrix, nullptr, "XY");
  HelperGeneral::CheckHistogramsForAxisIdentity(histoRespMatrix, histoYieldExt, "X");

  // ------------- expo parameters estimate --------------------------------------------
  TH1* histEff = histoRespMatrix->ProjectionX(); // for expo parameters estimate only
  TH1* histPreliminaryCorrectedYield = dynamic_cast<TH1*>(histoYieldExt->Clone());
  histPreliminaryCorrectedYield->SetDirectory(nullptr);
  histPreliminaryCorrectedYield->Divide(histEff);
  const double preliminaryYield = histPreliminaryCorrectedYield->Integral();
  const double preliminaryA = preliminaryYield / HelperGeneral::LifetimeLambdaC / (std::exp(-lo / HelperGeneral::LifetimeLambdaC) - std::exp(-hi / HelperGeneral::LifetimeLambdaC));
  // -----------------------------------------------------------------------------------

  ForwardFoldedExpo ffe{};
  ffe.SetResponseMatrix(histoRespMatrix);
  ffe.Init();

  TGraphErrors grTau;
  TGraphErrors grA;
  TGraphErrors grChi2;
  TGraphErrors grEta;
  grTau.SetName("grTau");
  grA.SetName("grA");
  grChi2.SetName("grChi2");
  grEta.SetName("grEta");
  grTau.GetYaxis()->SetTitle("#tau (fs)");
  grA.GetYaxis()->SetTitle("A");
  grChi2.GetYaxis()->SetTitle("#chi^{2} / ndf");
  grEta.GetYaxis()->SetTitle("#eta, %");

  for(const auto& gr : {&grTau, &grA, &grChi2, &grEta}) {
    gr->GetXaxis()->SetTitle("drop set");
  }

  const int nBins = histoYield->GetNbinsX();
  const int nDropSets = RunMode > HelperCtFit::RunModes::MeanFitOnly ? 1 << nBins : 1;
  double referenceSigmaTau{HelperGeneral::UndefValueDouble};

  auto FitLifetime = [&](TH1* histoRec, int dropSet) {
    const auto binsToDrop = HelperCtFit::EvalBinsToDrop(dropSet);

    TF1* fitFunc = new TF1("fitFunc", ffe, lo, hi, 2);
    fitFunc->SetParameters(preliminaryA, HelperGeneral::LifetimeLambdaC);
    fitFunc->SetNpx(1000);

    histoRec->Fit(fitFunc, "", "", lo, hi);

    grA.AddPoint(dropSet, fitFunc->GetParameter(0)*1000);
    grA.SetPointError(grA.GetN()-1, 0., fitFunc->GetParError(0)*1000);
    grTau.AddPoint(dropSet, fitFunc->GetParameter(1)*1000);
    grTau.SetPointError(grTau.GetN()-1, 0., fitFunc->GetParError(1)*1000);

    const auto chi2 = fitFunc->GetChisquare();
    const auto ndf = fitFunc->GetNDF();
    if(ndf != 0) grChi2.AddPoint(dropSet, chi2 / ndf);

    if (dropSet == 0) referenceSigmaTau = fitFunc->GetParError(1);
    if(binsToDrop.size() < 2) {
      const double currentSigmaTau = fitFunc->GetParError(1);
      const double eta = (currentSigmaTau / referenceSigmaTau - 1.) * 100.;
      grEta.AddPoint(dropSet, eta);
    }

    const std::pair<std::string, std::string> fitResults = HelperCtFit::FitResults(fitFunc);

    constexpr float textX1 = 0.70;
    constexpr float textX2 = 0.86;
    constexpr float textY2 = 0.86;
    constexpr float textYStep = 0.06;

    TCanvas ccFit("ccFit", "");
    ccFit.SetCanvasSize(1200, 800);
    histoRec->Draw("E1 same");
    std::string droppedBins{};
    for(const auto& binToDrop : binsToDrop) {
      droppedBins.append(std::to_string(binToDrop));
      droppedBins.append(", ");
    }
    if(dropSet != 0) {
      droppedBins.pop_back();
      droppedBins.pop_back();
      histoRec->SetTitle(("drop set = " + std::to_string(dropSet) + ", dropped bins: {" + droppedBins + "}").c_str());
    }
    fitFunc->Draw("same");
    HelperPlot::AddOneLineText(fitResults.first, {textX1, textY2 - 3*textYStep, textX2, textY2 - 2*textYStep}, "brNDC", 0.04);
    HelperPlot::AddOneLineText(fitResults.second, {textX1, textY2 - 4*textYStep, textX2, textY2 - 3*textYStep}, "brNDC", 0.04);
    HelperPlot::AddOneLineText(HelperCtFit::LifetimePdgSv, {textX1, textY2 - 6*textYStep, textX2, textY2 - 5*textYStep}, "brNDC", 0.04);
    ccFit.Print("ctfit.pdf", "pdf");

    TH1* hRatio = dynamic_cast<TH1*>(histoRec->Clone());
    HelperPlot::ScalePlotVertically(hRatio, histoRec, 2);
    hRatio->GetYaxis()->SetTitle("Data / Fit");
    HelperMath::DivideHistoByFunction(hRatio, fitFunc, "I");
    hRatio->SetBinContent(1, hRatio->GetBinContent(2)); // To exclude the left-most bin from CustomizeHistogramsYRange
    HelperPlot::CustomizeHistogramsYRange({hRatio}, false, -1e9, 1e9, 0.7);
    hRatio->SetBinContent(1, 0.); // revert 2 lines above workaround
    TCanvas ccRatio("ccRatio", "");
    HelperPlot::ScaleCanvasVertically(&ccRatio, &ccFit, 2);
    TF1 oneline("oneline", "[0]", histoRec->GetXaxis()->GetXmin(), histoRec->GetXaxis()->GetXmax());
    oneline.SetParameter(0, 1);
    oneline.SetLineColor(kBlack);
    oneline.SetLineStyle(7);
    hRatio->Draw("HIST PE same");
    oneline.Draw("same");
    ccRatio.Print("ctfit.pdf", "pdf");

    TH1* hResidual = dynamic_cast<TH1*>(histoRec->Clone());
    HelperPlot::ScalePlotVertically(hResidual, histoRec, 2);
    hResidual->GetYaxis()->SetTitle("(Data - Fit) / #sigma_{Data}");
    HelperMath::EvalNormDifferenceHistoFromFunction(hResidual, fitFunc, "I");
    hResidual->SetBinContent(1, hResidual->GetBinContent(2)); // To exclude the left-most bin from CustomizeHistogramsYRange
    HelperPlot::CustomizeHistogramsYRange({hResidual}, false, -1e9, 1e9, 0.7);
    hResidual->SetBinContent(1, 0.); // revert 2 lines above workaround
    TCanvas ccResidual("ccResidual", "");
    HelperPlot::ScaleCanvasVertically(&ccResidual, &ccFit, 2);
    TF1 zeroline("oneline", "[0]", 0, 2);
    zeroline.SetParameter(0, 0);
    zeroline.SetLineColor(kBlack);
    zeroline.SetLineStyle(7);
    hResidual->Draw("HIST PE same");
    zeroline.Draw("same");
    ccResidual.Print("ctfit.pdf", "pdf");

    TFile* fileOut{};
    if constexpr (IsSaveCanvasAsRoot) {
      fileOut = TFile::Open("ctfit.root", "recreate");
      fitFunc->Write();
      ccFit.Write();
      ccRatio.Write();
      ccResidual.Write();
      fileOut->Close();
    }
  };

  TCanvas emptycanvas("emptycanvas", "", 1200, 800);
  emptycanvas.Print("ctfit.pdf[", "pdf");

  for(int iDropSet=0; iDropSet<nDropSets; ++iDropSet) {
    const auto binsToDrop = HelperCtFit::EvalBinsToDrop(iDropSet);
    const auto binsToDropSize = static_cast<int>(binsToDrop.size());
    if(RunMode == HelperCtFit::RunModes::AllPossible && binsToDropSize >= nBins-1) continue;
    if(RunMode == HelperCtFit::RunModes::AllWoOne && binsToDropSize > 1) continue;
    TH1* hReco = dynamic_cast<TH1*>(histoYieldExt->Clone());

    for(const auto& binToDrop : binsToDrop) {
      HelperCtFit::ExcludeBin(hReco, binToDrop+1); // +1 because of additionalBinEdges
    }

    FitLifetime(hReco, iDropSet);
  } // nDropSets

  emptycanvas.Print("ctfit.pdf]", "pdf");

  HelperCtFit::PrintCanvasDrops(&grTau, nBins, "(");
  HelperCtFit::PrintCanvasDrops(&grA, nBins, "");
  HelperCtFit::PrintCanvasDrops(&grChi2, nBins, "");
  HelperCtFit::PrintCanvasDrops(&grEta, nBins, ")");

  fileRespMatrix->Close();
  fileYield->Close();
}

TH1* ExtendHistoWithEmptyBinsLeft(TH1* histo, const std::vector<double>& additionalBinEdges) {
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
  return hResult;
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
