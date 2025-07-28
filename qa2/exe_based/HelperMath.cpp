//
// Created by oleksii on 21.07.25.
//

#include "HelperMath.hpp"

#include "HelperGeneral.hpp"

#include <TF1.h>
#include <TH1.h>

#include <stdexcept>

std::pair<float, float> HelperMath::EstimateExpoParameters(TH1* h, float lo, float hi) {
  const int ilo = h->FindBin(lo);
  const int ihi = h->FindBin(hi);
  const float flo = h->GetBinContent(ilo)/* * h->GetBinWidth(ilo)*/;
  const float fhi = h->GetBinContent(ihi)/* * h->GetBinWidth(ihi)*/;
  const float tau = (hi-lo)/std::log(flo/fhi);
  const float A = flo / std::exp(-lo/tau);
  return std::make_pair(A, tau);
}

std::pair<double, double> HelperMath::DetermineWorkingRangesTH1(const TH1* histo, double leftMargin, double rightMargin) {
  double left{-999.};
  double right = {-999.};
  double leftTail = 0.;
  double rightTail = 0.;

  const double integral = histo->Integral(0, histo->GetNbinsX()+1);
  const int nBins = histo->GetNbinsX();

  for(int iBin=0; iBin<=nBins+1; iBin++) {
    leftTail += histo->GetBinContent(iBin);
    if(leftTail > integral*leftMargin) {
      left = iBin == 0 ? -std::numeric_limits<double>::infinity() : histo->GetBinLowEdge(iBin);
      break;
    }
  }

  for(int iBin=nBins+1; iBin>=0; iBin--) {
    rightTail += histo->GetBinContent(iBin);
    if(rightTail > integral*rightMargin) {
      right = iBin == nBins+1 ? std::numeric_limits<double>::infinity() : histo->GetBinLowEdge(iBin+1);
      break;
    }
  }

  return std::make_pair(left, right);
}

TGraph* HelperMath::EvaluateMovingAverage(const TGraph* graphIn, int aveLength, bool excludeOwnPoint) {
  TGraph* graphOut = new TGraph();
  EvaluateMovingAverage(graphIn, graphOut, aveLength);

  return graphOut;
}

void HelperMath::EvaluateMovingAverage(const TGraph* graphIn, TGraph* graphOut, int radius, bool isExcludeOwnPoint) {
  if(graphIn == nullptr || graphOut == nullptr) throw std::runtime_error("HelperMath::EvaluateMovingAverage(): graphIn == nullptr || graphOut == nullptr");

  const int nPoints = graphIn->GetN();
  for(int iPoint=0; iPoint<nPoints; ++iPoint) {
    int nLocalPoints{0};
    double value{0.};
    const int localRadius = std::min(radius, std::min(iPoint, nPoints-1-iPoint));
    for(int iLocalPoint=iPoint-localRadius; iLocalPoint<=iPoint+localRadius; ++iLocalPoint) {
      if(isExcludeOwnPoint && iLocalPoint == iPoint && localRadius!=0) continue;
      value += graphIn->GetPointY(iLocalPoint);
      ++nLocalPoints;
    }
    graphOut->SetPoint(iPoint, graphIn->GetPointX(iPoint), value/nLocalPoints);
  } // nPoints
}

void HelperMath::DivideGraph(TGraph* num, const TGraph* den) {
  if(num == nullptr || den == nullptr) throw std::runtime_error("HelperMath::DivideGraph(): num == nullptr || den == nullptr");
  const int nPoints = num->GetN();
  if(den->GetN() != nPoints) throw std::runtime_error("HelperMath::DivideGraph(): den->GetN() != nPoints");
  for(int iPoint=0; iPoint<nPoints; ++iPoint) {
    const double xNum = num->GetPointX(iPoint);
    const double xDen = den->GetPointX(iPoint);
    const double yNum = num->GetPointY(iPoint);
    const double yDen = den->GetPointY(iPoint);
    if(std::fabs(xNum - xDen)>1e-4) throw std::runtime_error("HelperMath::DivideGraph(): num and den points X coordinates do not match");
    num->SetPointY(iPoint, yNum/yDen);
  }
}

TF1* HelperMath::FitLifetimeHisto(TH1* histo, const std::string& option) {
  const double lo = histo->GetBinLowEdge(1) + 1e-3;
  const double hi = histo->GetBinLowEdge(histo->GetNbinsX()+1) - 1e-3;
  auto parEst = HelperMath::EstimateExpoParameters(histo, lo, hi);
  TF1* fitFunc = new TF1("fitFunc", "[0]*TMath::Exp(-x/[1])", lo, hi);
  fitFunc->SetParameters(parEst.first, parEst.second);
  histo->Fit(fitFunc, ("0"+option).c_str(), "", lo, hi);
  fitFunc->SetLineColor(histo->GetLineColor());

  return fitFunc;
}

void HelperMath::DivideFunctionByHisto(TH1* histo, TF1* func, const std::string& option) {
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
    throw std::runtime_error("HelperMath::DivideFunctionByHisto() - 'option' must be either empty string or I");
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

TH1* HelperMath::MergeHistograms(const std::vector<TH1*>& histos) {
  for(const auto& h : histos) {
    HelperGeneral::CheckHistogramsForXaxisIdentity(h, histos.at(0));
  }

  TH1* hResult = dynamic_cast<TH1*>(histos.at(0)->Clone("merged"));
  hResult->Sumw2();
  hResult->SetDirectory(nullptr);
  for(size_t iH=1, nHs=histos.size(); iH<nHs; ++iH) {
    hResult->Add(histos.at(iH));
  }

  return hResult;
}
