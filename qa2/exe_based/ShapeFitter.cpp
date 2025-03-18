//
// Created by oleksii on 04.02.25.
//

#include "ShapeFitter.hpp"

#include "Helper.hpp"
#include "Shapes.hpp"

using namespace Helper;

void ShapeFitter::SetSideBands(double le, double li, double ri, double re) {
  left_sideband_external_ = le;
  left_sideband_internal_ = li;
  right_sideband_internal_ = ri;
  right_sideband_external_ = re;
}

void ShapeFitter::SetSideBandsViaSigma(double le, double li, double ri, double re) {
  if(expected_mu_ == UndefValueDouble || expected_sigma_ == UndefValueDouble) {
    throw std::runtime_error("SetSideBandsViaSigma() - expected_mu_ and expected_sigma_ must be set first");
  }
  left_sideband_external_ = expected_mu_ + le*expected_sigma_;
  left_sideband_internal_ = expected_mu_ + li*expected_sigma_;
  right_sideband_internal_ = expected_mu_ + ri*expected_sigma_;
  right_sideband_external_ = expected_mu_ + re*expected_sigma_;
}

void ShapeFitter::Fit() {
  PrepareHistoSidebands();
  FitSideBands(histo_sidebands_);
  PrepareHistoPeak();
  FitPeak(histo_peak_, peak_shape_);
  DefineAll(left_sideband_external_, right_sideband_external_);
}

void ShapeFitter::FitPeak(TH1D* h, const std::string& peakFunc) {
  DefinePeak(h, left_sideband_external_, right_sideband_external_);
  peak_fit_->SetNpx(1000);

  h->Fit(peak_fit_, "R0");
  chi2_peak_ = peak_fit_->GetChisquare();
  ndf_peak_ = peak_fit_->GetNDF();
}

void ShapeFitter::FitSideBands(TH1D* h) {
  DefineSideBand(left_sideband_external_, right_sideband_external_);
  sidebands_fit_->SetNpx(1000);

  h->Fit(sidebands_fit_, "R0");
  chi2_sideband_ = sidebands_fit_->GetChisquare();
  ndf_sideband_ = sidebands_fit_->GetNDF();
}

TPaveText* ShapeFitter::ConvertFitParametersToText(const std::string& funcType, std::array<float, 2> coordinatesLeftUpperCorner) const {
  TF1* func{nullptr};
  double chi2;
  int ndf;
  if(funcType == "peak") {
    func = peak_fit_;
    chi2 = chi2_peak_;
    ndf = ndf_peak_;
  } else if(funcType == "bg") {
    func = sidebands_fit_;
    chi2 = chi2_sideband_;
    ndf = ndf_sideband_;
  } else if(funcType == "all") {
    func = all_fit_;
    chi2 = chi2_all_;
    ndf = ndf_all_;
  }
  else throw std::runtime_error("ShapeFitter::ConvertFitParametersToText() - funcType must be either peak or bg");
  const float x1 = coordinatesLeftUpperCorner.at(0);
  const float y2 = coordinatesLeftUpperCorner.at(1);
  const float x2 = x1 + 0.2;
  const float y1 = y2 - 0.05 * func->GetNpar();
  auto* ptpar = new TPaveText(x1, y1, x2, y2, "brNDC");
  ptpar->SetFillColor(0);
  ptpar->SetTextSize(0.025);
  ptpar->SetTextFont(22);
  ptpar->AddText(func->GetTitle());
  ptpar->AddText(("#chi^{2}/ndf (" + funcType + ") = " + to_string_with_precision(chi2, 2) + " / " + std::to_string(ndf)).c_str());
  const int nPar = func->GetNpar();
  for(int iPar=0; iPar<nPar; iPar++) {
    ptpar->AddText((static_cast<std::string>(func->GetParName(iPar)) + " = " +
                  to_string_with_significant_figures(func->GetParameter(iPar), 3) + " #pm " +
                  to_string_with_significant_figures(func->GetParError(iPar), 3)).c_str());
  }

  return ptpar;
}

void ShapeFitter::PrepareHistoSidebands() {
  histo_sidebands_ = dynamic_cast<TH1D*>(histo_in_->Clone());
  const int nBins = histo_sidebands_->GetNbinsX();
  for(int iBin = 1; iBin <= nBins; iBin++) {
    const double binCenter = histo_sidebands_->GetBinCenter(iBin);
    if(binCenter < left_sideband_external_ ||
       (binCenter > left_sideband_internal_ && binCenter < right_sideband_internal_) ||
       binCenter > right_sideband_external_) {
      histo_sidebands_->SetBinContent(iBin, 0.);
      histo_sidebands_->SetBinError(iBin, 0.);
    }
  }
}

void ShapeFitter::PrepareHistoPeak() {
  histo_peak_ = dynamic_cast<TH1D*>(histo_in_->Clone());
  histo_peak_->Sumw2();
  histo_peak_->Add(sidebands_fit_, -1);

  const int nBins = histo_peak_->GetNbinsX();
  for (int iBin = 1; iBin <= nBins; iBin++) {
    const double binCenter = histo_peak_->GetBinCenter(iBin);
    if (binCenter < left_sideband_external_ || binCenter > right_sideband_external_) {
      histo_peak_->SetBinContent(iBin, 0.);
      histo_peak_->SetBinError(iBin, 0.);
    }
  }
}

void ShapeFitter::DefineSideBand(double left, double right) {
  if(bg_pol_n_ < -1 || bg_pol_n_ > PolN::nPars-2) {
    throw std::runtime_error("ShapeFitter::DefineSideBand() - bg_pol_n_ must be from 0 to " + std::to_string(PolN::nPars-2));
  }
  sidebands_fit_ = new TF1("sideband_fit", PolN::Shape, left, right, PolN::nPars);
  sidebands_fit_->FixParameter(PolN::kShift, expected_mu_);
  for(int iPar = bg_pol_n_+2; iPar < PolN::nPars; iPar++) {
    sidebands_fit_->FixParameter(iPar, 0.);
  }
  sidebands_fit_->SetParNames("#mu_{ref}", "a_{0}", "a_{1}", "a_{2}", "a_{3}");
  sidebands_fit_->SetTitle(("Pol" + std::to_string(bg_pol_n_)).c_str());
}

void ShapeFitter::DefineAll(double left, double right) {
  int nParsPeak{UndefValueInt};
  if(peak_shape_ == "Gaus") nParsPeak = Gaus::nPars;
  if(peak_shape_ == "DoubleGaus") nParsPeak = DoubleGaus::nPars;
  if(peak_shape_ == "DSCB") nParsPeak = DoubleSidedCrystalBall::nPars;
  auto AllShape = [&] (const double* x, const double* par) {
    if(peak_shape_ == "Gaus") return PolN::Shape(x, par) + Gaus::Shape(x, &par[PolN::nPars]);
    if(peak_shape_ == "DoubleGaus") return PolN::Shape(x, par) + DoubleGaus::Shape(x, &par[PolN::nPars]);
    if(peak_shape_ == "DSCB") return PolN::Shape(x, par) + DoubleSidedCrystalBall::Shape(x, &par[PolN::nPars]);
    throw std::runtime_error("ShapeFitter::DefineAll() - peak_shape_ must be one of the availables");
  };
  std::cout << "nParsPeak = " << nParsPeak << "\n";
  all_fit_ = new TF1("all_fit", AllShape, left, right, nParsPeak + PolN::nPars);
  CopyPasteParametersToAll(PolN::nPars, nParsPeak);
}

void ShapeFitter::CopyPasteParametersToAll(int nParsSideBand, int nParsPeak) {
  for (int iPar = 0; iPar < nParsSideBand; iPar++) {
    all_fit_->SetParameter(iPar, sidebands_fit_->GetParameter(iPar));
    all_fit_->SetParName(iPar, sidebands_fit_->GetParName(iPar));
    double minlimit, maxlimit;
    sidebands_fit_->GetParLimits(iPar, minlimit, maxlimit);
    all_fit_->SetParLimits(iPar, minlimit, maxlimit);
  }
  for (int iPar = 0; iPar < nParsPeak; iPar++) {
    all_fit_->SetParameter(nParsSideBand + iPar, peak_fit_->GetParameter(iPar));
    all_fit_->SetParName(nParsSideBand + iPar, peak_fit_->GetParName(iPar));
    double minlimit, maxlimit;
    peak_fit_->GetParLimits(iPar, minlimit, maxlimit);
    all_fit_->SetParLimits(nParsSideBand + iPar, minlimit, maxlimit);
  }
}

void ShapeFitter::DefinePeak(TH1D* histo, float left, float right) {
  if      (peak_shape_ == "Gaus")       DefinePeakGaus(histo, left, right);
  else if (peak_shape_ == "DSCB")       DefinePeakDSCB(histo, left, right);
  else if (peak_shape_ == "DoubleGaus") DefinePeakDoubleGaus(histo, left, right);
  else                                  throw std::runtime_error("ShapeFitter::DefinePeak(): peak_shape_ must be one of the available");
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

