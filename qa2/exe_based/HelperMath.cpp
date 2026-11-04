//
// Created by oleksii on 21.07.25.
//

#include "HelperMath.hpp"

#include "HelperGeneral.hpp"

#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TMatrixD.h>

#include <stdexcept>

using namespace HelperGeneral;

std::pair<double, double> HelperMath::EstimateExpoParameters(TH1* h) {
  int ilo{1};
  while(h->GetBinContent(ilo) == 0.) {
    ++ilo;
  }
  int ihi{h->GetNbinsX()};
  while(h->GetBinContent(ihi) == 0.) {
    --ihi;
  }
  const double lo = h->GetBinCenter(ilo);
  const double hi = h->GetBinCenter(ihi);
  const double flo = h->GetBinContent(ilo);
  const double fhi = h->GetBinContent(ihi);
  const double tau = (hi-lo)/std::log(flo/fhi);
  const double A = flo / std::exp(-lo/tau);
  return std::make_pair(A, tau);
}

TF1* HelperMath::FitLifetimeHisto(TH1* histo, const std::string& option) {
  const double lo = histo->GetBinLowEdge(1) + 1e-3;
  const double hi = histo->GetBinLowEdge(histo->GetNbinsX()+1) - 1e-3;
  const auto parEst = HelperMath::EstimateExpoParameters(histo);
  TF1* fitFunc = new TF1("fitFunc", "[0]*TMath::Exp(-x/[1])", lo, hi);
  fitFunc->SetParameters(parEst.first, parEst.second);
  histo->Fit(fitFunc, ("0"+option).c_str(), "", lo, hi);
  fitFunc->SetLineColor(histo->GetLineColor());

  return fitFunc;
}

void HelperMath::DivideHistoByFunction(TH1* histo, TF1* func, const std::string& option) {
  if(option.empty()) {
    histo->Divide(func);
  } else if(option == "I") {
    for(int iBin=1, nBins=histo->GetNbinsX(); iBin<=nBins; iBin++) {
      const double histoValue = histo->GetBinContent(iBin);
      const double histoError = histo->GetBinError(iBin);
      const double lo = histo->GetBinLowEdge(iBin);
      const double hi = histo->GetBinLowEdge(iBin+1);
      const double funcAverage = func->Integral(lo, hi) / (hi-lo);
      histo->SetBinContent(iBin, histoValue/funcAverage);
      histo->SetBinError(iBin, histoError/funcAverage);
    }
  } else {
    throw std::runtime_error("HelperMath::DivideHistoByFunction() - 'option' must be either empty string or I");
  }
}

void HelperMath::EvalNormDifferenceHistoFromFunction(TH1* histo, TF1* func, const std::string& option) {
  const bool isIntegral = option == "I" ? true : option.empty() ? false : throw std::runtime_error("HelperMath::EvalNormDifferenceHistoFromFunction() - 'option' must be either empty string or I");
  for(int iBin=1, nBins=histo->GetNbinsX(); iBin<=nBins; iBin++) {
    const double histoValue = histo->GetBinContent(iBin);
    const double histoError = histo->GetBinError(iBin);
    const double lo = histo->GetBinLowEdge(iBin);
    const double hi = histo->GetBinLowEdge(iBin+1);
    const double ce = histo->GetBinCenter(iBin);
    const double funcValue = isIntegral ? func->Integral(lo, hi) / (hi-lo) : func->Eval(ce);
    histo->SetBinContent(iBin, histoError != 0. ? (histoValue - funcValue) / histoError : 0.);
    histo->SetBinError(iBin, histoError != 0. ? 1e-9 : 0.);
  }
}

void HelperMath::InvertHisto(TH1* histo) {
  Sumw2IfNotYet(histo);
  for(int iBin=1, nBins=histo->GetNbinsX(); iBin<nBins; ++iBin) {
    const double value = histo->GetBinContent(iBin);
    const double error = histo->GetBinError(iBin);
    histo->SetBinContent(iBin, 1./value);
    histo->SetBinError(iBin, error / value / value);
  }
}

std::pair<TH1*, TH1*> HelperMath::EvaluateEfficiencyHisto(TH1* hNum, TH1* hDen) {
  const int nBins = hNum->GetNbinsX();
  HelperGeneral::CheckHistogramsForXaxisIdentity(hNum, hDen);

  TH1* hEff = dynamic_cast<TH1*>(hNum->Clone());
  TH1* hRelErr = dynamic_cast<TH1*>(hNum->Clone());
  hEff->Reset();
  hRelErr->Reset();

  auto EvalEfficiency = [](double num, double den) {
      if (den == 0.) return 0.;
      else return num / den;
  };

  auto EvalRelErrOfEfficiency = [](double num, double den) {
      if (num == 0. || den == 0.) return 0.;
      if (num > den) return 1.;

      return std::sqrt(1. / num - 1. / den);
  };

  auto EvalAbsErrOfRelErrOfEfficiency = [&](double num, double den) {
      if (num == 0. || den == 0. || num > den) return 0.;

      auto relErr = EvalRelErrOfEfficiency(num, den);
      return 1. / 2. / relErr * std::sqrt(1. / num / num / num + 1. / den / den / den - 2. / num / den / den);
  };

  for (int iBin = 1; iBin <= nBins; iBin++) {
    const double num = hNum->GetBinContent(iBin);
    const double den = hDen->GetBinContent(iBin);
    const double eff = EvalEfficiency(num, den);
    const double relErr = EvalRelErrOfEfficiency(num, den);
    const double absErrOnRelErr = EvalAbsErrOfRelErrOfEfficiency(num, den);
    hEff->SetBinContent(iBin, eff);
    hEff->SetBinError(iBin, eff * relErr);
    hRelErr->SetBinContent(iBin, relErr);
    hRelErr->SetBinError(iBin, absErrOnRelErr);
  }

  hEff->GetYaxis()->SetTitle("#varepsilon");
  hEff->SetTitle("");

  hRelErr->GetYaxis()->SetTitle("#varepsilon_{#varepsilon}");
  hRelErr->SetTitle("");

  return std::make_pair(hEff, hRelErr);
}

std::pair<TH2*, TH2*> HelperMath::EvaluateResponseMatrix(TH2* hRec, TH1* hGen) {
  HelperGeneral::CheckHistogramsForAxisIdentity<TH2, TH2>(hRec, nullptr, "XY");
  HelperGeneral::CheckHistogramsForXaxisIdentity(hRec, hGen);

  TH2* hResp = dynamic_cast<TH2*>(hRec->Clone());
  TH2* hRelErr = dynamic_cast<TH2*>(hRec->Clone());
  hResp->Reset();
  hRelErr->Reset();

  auto EvalResponce = [](double num, double den) {
    if (den == 0.) return 0.;
    else return num / den;
  };

  auto EvalRelErrOfResponce = [](double num, double den) {
    if (num == 0. || den == 0.) return 0.;
    if (num > den) {
      std::cout << "Warning! HelperMath::EvaluateResponseMatrix::EvalRelErrOfResponce() num > den (" << num << " vs " << den << ")\n";
      return 1.;
    }

    return std::sqrt(1. / num - 1. / den);
  };

  auto EvalAbsErrOfRelErrOfResponce = [&](double num, double den) {
    if (num == 0. || den == 0. || num > den) return 0.;

    auto relErr = EvalRelErrOfResponce(num, den);
    return 1. / 2. / relErr * std::sqrt(1. / num / num / num + 1. / den / den / den - 2. / num / den / den);
  };

  const int nBins = hRec->GetNbinsX();
  for(int jBin=1; jBin<= nBins; ++jBin) {
    const double den = hGen->GetBinContent(jBin);
    for(int iBin=1; iBin<=nBins; ++iBin) {
      const double num = hRec->GetBinContent(iBin, jBin);
      const double resp = EvalResponce(num, den);
      const double relErr = EvalRelErrOfResponce(num, den);
      const double absErrOnRelErr = EvalAbsErrOfRelErrOfResponce(num, den);
      hResp->SetBinContent(iBin, jBin, resp);
      hResp->SetBinError(iBin, jBin, resp * relErr);
      hRelErr->SetBinContent(iBin, jBin, relErr);
      hRelErr->SetBinError(iBin, jBin, absErrOnRelErr);
    } // iBin, rec
  } // jBin, gen

  return std::make_pair(hResp, hRelErr);
}

double HelperMath::EvalErrorFitFunction(double x, TF1* func, const TMatrixDSym& cov) {
  const int nPars = func->GetNpar();
  TMatrixD dfdp(nPars, 1);
  for (int iPar = 0; iPar < nPars; iPar++) {
    dfdp[iPar][0] = func->GradientPar(iPar, &x);
  }
  TMatrixD dfdp_T = dfdp;
  dfdp_T.T();

  double result = std::sqrt((dfdp_T * cov * dfdp)[0][0]);
  if(!std::isfinite(result)) {
    result = 0.;
  }

  return result;
}

TH1* HelperMath::CutSubHistogram(const TH1* histoIn, double lo, double hi) {
  if(lo >= hi) throw std::runtime_error("HelperMath::CutSubHistogram(): lo >= hi");

  const double tolerance = 1e-6;
  int binLoIn{UndefValueInt};
  bool isEndReached{false};
  std::vector<double> binEdges;
  for(int iBin=1, nBins=histoIn->GetNbinsX(); iBin<=nBins+1; ++iBin) {
    const double binLowEdge = histoIn->GetBinLowEdge(iBin);
    if(EqualFloating(binLowEdge, lo, tolerance)) binLoIn = iBin;
    if(binLoIn != UndefValueInt) binEdges.emplace_back(binLowEdge);
    if(EqualFloating(binLowEdge, hi, tolerance)) {
      isEndReached = true;
      break;
    }
  } // histoIn bins
  if(binLoIn == UndefValueInt || !isEndReached) throw std::runtime_error("HelperMath::CutSubHistogram(): either lo or hi does not match any of histoIn bin edges");

  TH1* histoOut = new TH1D("", "", binEdges.size()-1, binEdges.data());
  histoOut->SetDirectory(nullptr);
  if(histoIn->GetSumw2N() > 0) histoOut->Sumw2();
  histoOut->GetXaxis()->SetTitle(histoIn->GetXaxis()->GetTitle());
  histoOut->GetYaxis()->SetTitle(histoIn->GetYaxis()->GetTitle());
  histoOut->SetName(histoIn->GetName());
  histoOut->SetTitle(histoIn->GetTitle());
  histoOut->SetLineColor(histoIn->GetLineColor());
  histoOut->SetLineWidth(histoIn->GetLineWidth());
  histoOut->SetLineStyle(histoIn->GetLineStyle());
  histoOut->SetMarkerColor(histoIn->GetMarkerColor());
  histoOut->SetMarkerStyle(histoIn->GetMarkerStyle());
  histoOut->SetMarkerSize(histoIn->GetMarkerSize());
  histoOut->SetFillColor(histoIn->GetFillColor());
  histoOut->SetFillStyle(histoIn->GetFillStyle());
  histoOut->SetOption(histoIn->GetOption());
  for(int iBin=1, nBins=binEdges.size()-1; iBin<=nBins; ++iBin) {
    const double value = histoIn->GetBinContent(binLoIn-1 + iBin);
    const double error = histoIn->GetBinError(binLoIn-1 + iBin);
    histoOut->SetBinContent(iBin, value);
    histoOut->SetBinError(iBin, error);
  }

  return histoOut;
}
