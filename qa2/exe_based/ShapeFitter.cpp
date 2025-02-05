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

  std::cout << "peak_fit_->GetName() = " << peak_fit_->GetName() << "\n";
  std::cout << "peak_fit_->GetTitle() = " << peak_fit_->GetTitle() << "\n";

  h->Fit(peak_fit_, "R");
  chi2_peak_ = peak_fit_->GetChisquare() / peak_fit_->GetNDF();
}

TPaveText ShapeFitter::FitParametersToText(float x1, float y1, float x2, float y2) const {
  TPaveText ptpar(x1, y1, x2, y2, "brNDC");
  ptpar.SetFillColor(0);
  ptpar.SetTextSize(0.025);
  ptpar.SetTextFont(22);
  ptpar.AddText(("#chi^{2}/_{ndf} (peak) = " + to_string_with_precision(chi2_peak_, 2)).c_str());
  const int nPar = peak_fit_->GetNpar();
  for(int iPar=0; iPar<nPar; iPar++) {
    ptpar.AddText((static_cast<std::string>(peak_fit_->GetParName(iPar)) + " = " +
                  to_string_with_precision(peak_fit_->GetParameter(iPar), 2) + " #pm " +
                  to_string_with_precision(peak_fit_->GetParError(iPar), 2)).c_str());
  }
}

void ShapeFitter::DefinePeakGaus(TH1D* histo, double left, double right) {
  const int Npar = 4;
  peak_fit_ = new TF1("peak_fit", Gaus, left, right, Npar);
  peak_fit_->SetParameter(0, histo->Interpolate(expected_mu_));
  peak_fit_->FixParameter(1, expected_mu_);
  peak_fit_->SetParameter(2, 0);
  peak_fit_->SetParameter(3, expected_sigma_);
  peak_fit_->SetParLimits(0, 0, 10 * histo->Interpolate(expected_mu_));
  peak_fit_->SetParLimits(2, -5*expected_sigma_, 5*expected_sigma_);
  peak_fit_->SetParLimits(3, 0, 5*expected_sigma_);
  peak_fit_->SetParNames("Height", "#mu_{ref}", "#mu - #mu_{ref}", "#sigma");
}

void ShapeFitter::DefinePeakDoubleGaus(TH1D* histo, double left, double right) {
  const int Npar = 6;
  peak_fit_ = new TF1("peak_fit", DoubleGaus, left, right, Npar);
  peak_fit_->SetParameter(0, histo->Interpolate(expected_mu_) / 2);
  peak_fit_->SetParameter(1, histo->Interpolate(expected_mu_) / 2);
  peak_fit_->FixParameter(2, expected_mu_);
  peak_fit_->SetParameter(3, 0);
  peak_fit_->SetParameter(4, 0.5*expected_sigma_);
  peak_fit_->SetParameter(5, 2*expected_sigma_);
  peak_fit_->SetParLimits(0, 0, 10 * histo->Interpolate(expected_mu_));
  peak_fit_->SetParLimits(1, 0, 10 * histo->Interpolate(expected_mu_));
  peak_fit_->SetParLimits(3, -5*expected_sigma_, 5*expected_sigma_);
  peak_fit_->SetParLimits(4, 0, 5*expected_sigma_);
  peak_fit_->SetParLimits(5, 0, 5*expected_sigma_);
  peak_fit_->SetParNames("Height1", "Height2", "#mu_{ref}", "#mu - #mu_{ref}", "#sigma_{1}", "#sigma_{2}");
}

void ShapeFitter::DefinePeakDSCB(TH1D* histo, float left, float right) {
  const int Npar = 8;
  peak_fit_ = new TF1("sgnl_fit", DoubleSidedCrystalBall, left, right, Npar);
  peak_fit_->SetParameter(0, histo->Interpolate(expected_mu_));
  peak_fit_->FixParameter(1, expected_mu_);
  peak_fit_->SetParameter(2, 0);
  peak_fit_->SetParameter(3, expected_sigma_);
  peak_fit_->SetParameter(4, 1.);
  peak_fit_->SetParameter(5, 1.);
  peak_fit_->SetParameter(6, 1.);
  peak_fit_->SetParameter(7, 1.);
  peak_fit_->SetParLimits(0, 0, 10 * histo->Interpolate(expected_mu_));
  peak_fit_->SetParLimits(4, 0, 100);
  peak_fit_->SetParLimits(5, -5, 5);
  peak_fit_->SetParLimits(6, 0, 100);
  peak_fit_->SetParLimits(7, -5, 5);
  peak_fit_->SetParNames("Height1", "#mu_{ref}", "#mu - #mu_{ref}", "#sigma", "a_{1}", "log_{10}n_{1}", "a_{2}", "log_{10}n_{2}");
}

