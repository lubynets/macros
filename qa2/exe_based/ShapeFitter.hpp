//
// Created by oleksii on 04.02.25.
//

#ifndef QA2_SHAPEFITTER_HPP
#define QA2_SHAPEFITTER_HPP

#include <TF1.h>
#include <TH1.h>
#include <TPaveText.h>

class ShapeFitter {
 public:
  explicit ShapeFitter(TH1D* histo) { histo_in_ = histo; }
  virtual ~ShapeFitter() = default;

  const TF1* const GetPeakFunc() const { return peak_fit_; }
  double GetPeakChi2() const { return chi2_peak_; }

  void SetExpectedMu(float value) { expected_mu_ = value; }
  void SetExpectedSigma(float value) { expected_sigma_ = value; }

  void Fit(const std::string& peakFunc = "Gaus");

  TPaveText FitParametersToText(float x1, float y1, float x2, float y2) const;

 private:
  void FitPeak(TH1D* h, const std::string& peakFunc);

  void DefinePeakGaus(TH1D* h, double left, double right);
  void DefinePeakDoubleGaus(TH1D* histo, double left, double right);
  void DefinePeakDSCB(TH1D* histo, float left, float right);

  TH1D* histo_in_{nullptr};
  TF1* peak_fit_{nullptr};
  double chi2_peak_{-999.};

  double expected_mu_{-999.};
  double expected_sigma_{-999.};
};
#endif //QA2_SHAPEFITTER_HPP
