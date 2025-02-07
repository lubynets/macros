//
// Created by oleksii on 04.02.25.
//

#include "ShapeFitter.hpp"

#include "Helper.hpp"
#include "Shapes.hpp"

using namespace Helper;

void ShapeFitter::Fit(const std::string& peakFunc) {
  TH1D* histoPeak = histo_in_; // TODO be fixed when BG shape is defined
  FitPeak(histoPeak, peakFunc);
}

void ShapeFitter::FitPeak(TH1D* h, const std::string& peakFunc) {
  const double leftRange = h->GetBinCenter(1);
  const double rightRange = h->GetBinCenter(h->GetNbinsX());
  if      (peakFunc == "Gaus")       DefinePeakGaus(h, leftRange, rightRange);
  else if (peakFunc == "DSCB")       DefinePeakDSCB(h, leftRange, rightRange);
  else if (peakFunc == "DoubleGaus") DefinePeakDoubleGaus(h, leftRange, rightRange);
  else                         throw std::runtime_error("ShapeFitter::FitPeak(): peakFunc must be one of the available");
  peak_fit_->SetNpx(1000);

  h->Fit(peak_fit_, "R");
  chi2_peak_ = peak_fit_->GetChisquare() / peak_fit_->GetNDF();
}

TPaveText* ShapeFitter::FitParametersToText(float x1, float y1, float x2, float y2) const {
  auto* ptpar = new TPaveText(x1, y1, x2, y2, "brNDC");
  ptpar->SetFillColor(0);
  ptpar->SetTextSize(0.025);
  ptpar->SetTextFont(22);
  ptpar->AddText(peak_fit_->GetTitle());
  ptpar->AddText(("#chi^{2}/ndf (peak) = " + to_string_with_precision(chi2_peak_, 2)).c_str());
  const int nPar = peak_fit_->GetNpar();
  for(int iPar=0; iPar<nPar; iPar++) {
    ptpar->AddText((static_cast<std::string>(peak_fit_->GetParName(iPar)) + " = " +
                  to_string_with_significant_figures(peak_fit_->GetParameter(iPar), 3) + " #pm " +
                  to_string_with_significant_figures(peak_fit_->GetParError(iPar), 3)).c_str());
  }

  return ptpar;
}

void ShapeFitter::DefinePeakGaus(TH1D* histo, double left, double right) {
  peak_fit_ = new TF1("peak_fit", Gaus::Shape, left, right, Gaus::nPars);
  peak_fit_->SetParameter(Gaus::kFactor, histo->Interpolate(expected_mu_));
  peak_fit_->FixParameter(Gaus::kShift, expected_mu_);
  peak_fit_->SetParameter(Gaus::kMu, 0);
  peak_fit_->SetParameter(Gaus::kSigma, expected_sigma_);
  peak_fit_->SetParLimits(Gaus::kFactor, 0, 10 * histo->Interpolate(expected_mu_));
  peak_fit_->SetParLimits(Gaus::kMu, -5*expected_sigma_, 5*expected_sigma_);
  peak_fit_->SetParLimits(Gaus::kSigma, 0, 5*expected_sigma_);
  peak_fit_->SetParNames("Height", "#mu_{ref}", "#mu - #mu_{ref}", "#sigma");
  peak_fit_->SetTitle("Gaus");
}

void ShapeFitter::DefinePeakDoubleGaus(TH1D* histo, double left, double right) {
  peak_fit_ = new TF1("peak_fit", DoubleGaus::Shape, left, right, DoubleGaus::nPars);
  peak_fit_->SetParameter(DoubleGaus::kFactor1, histo->Interpolate(expected_mu_) / 2);
  peak_fit_->SetParameter(DoubleGaus::kFactor2, histo->Interpolate(expected_mu_) / 2);
  peak_fit_->FixParameter(DoubleGaus::kShift, expected_mu_);
  peak_fit_->SetParameter(DoubleGaus::kMu, 0);
  peak_fit_->SetParameter(DoubleGaus::kSigma1, 0.5*expected_sigma_);
  peak_fit_->SetParameter(DoubleGaus::kSigma2, 2*expected_sigma_);
  peak_fit_->SetParLimits(DoubleGaus::kFactor1, 0, 10 * histo->Interpolate(expected_mu_));
  peak_fit_->SetParLimits(DoubleGaus::kFactor2, 0, 10 * histo->Interpolate(expected_mu_));
  peak_fit_->SetParLimits(DoubleGaus::kMu, -5*expected_sigma_, 5*expected_sigma_);
  peak_fit_->SetParLimits(DoubleGaus::kSigma1, 0, 5*expected_sigma_);
  peak_fit_->SetParLimits(DoubleGaus::kSigma2, 0, 5*expected_sigma_);
  peak_fit_->SetParNames("Height1", "Height2", "#mu_{ref}", "#mu - #mu_{ref}", "#sigma_{1}", "#sigma_{2}");
  peak_fit_->SetTitle("DoubleGaus");
}

void ShapeFitter::DefinePeakDSCB(TH1D* histo, float left, float right) {
  peak_fit_ = new TF1("sgnl_fit", DoubleSidedCrystalBall::Shape, left, right, DoubleSidedCrystalBall::nPars);
  peak_fit_->SetParameter(DoubleSidedCrystalBall::kFactor, histo->Interpolate(expected_mu_));
  peak_fit_->FixParameter(DoubleSidedCrystalBall::kShift, expected_mu_);
  peak_fit_->SetParameter(DoubleSidedCrystalBall::kMu, 0);
  peak_fit_->SetParameter(DoubleSidedCrystalBall::kSigma, expected_sigma_);
  peak_fit_->SetParameter(DoubleSidedCrystalBall::kA1, 1.);
  peak_fit_->SetParameter(DoubleSidedCrystalBall::kN1, 1.);
  peak_fit_->SetParameter(DoubleSidedCrystalBall::kA2, 1.);
  peak_fit_->SetParameter(DoubleSidedCrystalBall::kN2, 1.);
  peak_fit_->SetParLimits(DoubleSidedCrystalBall::kFactor, 0, 10 * histo->Interpolate(expected_mu_));
  peak_fit_->SetParLimits(DoubleSidedCrystalBall::kMu, -5*expected_sigma_, 5*expected_sigma_);
  peak_fit_->SetParLimits(DoubleSidedCrystalBall::kSigma, 0, 5*expected_sigma_);
  peak_fit_->SetParLimits(DoubleSidedCrystalBall::kA1, 0, 100);
  peak_fit_->SetParLimits(DoubleSidedCrystalBall::kN1, -5, 5);
  peak_fit_->SetParLimits(DoubleSidedCrystalBall::kA2, 0, 100);
  peak_fit_->SetParLimits(DoubleSidedCrystalBall::kN2, -5, 5);
  peak_fit_->SetParNames("Height", "#mu_{ref}", "#mu - #mu_{ref}", "#sigma", "a_{1}", "log_{10}n_{1}", "a_{2}", "log_{10}n_{2}");
  peak_fit_->SetTitle("DSCB");
}

