//
// Created by oleksii on 04.02.25.
//

#ifndef QA2_SHAPEFITTER_HPP
#define QA2_SHAPEFITTER_HPP

#include "Helper.hpp"

#include <TF1.h>
#include <TH1.h>
#include <TPaveText.h>

class ShapeFitter {
 public:
  explicit ShapeFitter(TH1D* histo) { histo_in_ = histo; }
  virtual ~ShapeFitter() = default;

  TF1* const GetPeakFunc() const { return peak_fit_; }
  double GetPeakChi2() const { return chi2_peak_; }
  TH1D* GetPeakHisto() const { return histo_peak_; }

  TF1* const GetSideBandFunc() const { return sidebands_fit_; }
  double GetSideBandChi2() const { return chi2_sideband_; }
  TH1D* GetSideBandHisto() const { return histo_sidebands_; }

  TF1* GetAllFunc() const { return all_fit_; }
  double GetAllChi2() const { return chi2_all_; }
  TH1D* GetAllHisto() const { return histo_in_; }

  void SetExpectedMu(double value) { expected_mu_ = value; }
  void SetExpectedSigma(double value) { expected_sigma_ = value; }
  void SetSideBands(double le, double li, double ri, double re);
  void SetSideBandsViaSigma(double le, double li, double ri, double re);
  void SetPeakShape(const std::string& shapeName) { peak_shape_ = shapeName; }
  void SetBgPolN(int polN) { bg_pol_n_ = polN; }

  void Fit();

  TPaveText* ConvertFitParametersToText(const std::string& funcType, std::array<float, 2> coordinatesLeftUpperCorner) const;

 private:
  void FitPeak(TH1D* h, const std::string& peakFunc);
  void FitSideBands(TH1D* h);

  void DefinePeakGaus(TH1D* h, double left, double right);
  void DefinePeakDoubleGaus(TH1D* histo, double left, double right);
  void DefinePeakDSCB(TH1D* histo, float left, float right);

  void CopyPasteParametersToAll(int nParsSideBand, int nParsPeak);

  void DefinePeak(TH1D* histo, float left, float right);
  void DefineSideBand(double left, double right);
  void DefineAll(double left, double right);

  void PrepareHistoSidebands();
  void PrepareHistoPeak();

  TH1D* histo_in_{nullptr};
  TH1D* histo_sidebands_{nullptr};
  TH1D* histo_peak_{nullptr};
  TF1* peak_fit_{nullptr};
  TF1* sidebands_fit_{nullptr};
  TF1* all_fit_{nullptr};

  double chi2_peak_{Helper::UndefValueDouble};
  double chi2_sideband_{Helper::UndefValueDouble};
  double chi2_all_{Helper::UndefValueDouble};
  int ndf_peak_{Helper::UndefValueInt};
  int ndf_sideband_{Helper::UndefValueInt};
  int ndf_all_{Helper::UndefValueInt};

  double expected_mu_{Helper::UndefValueDouble};
  double expected_sigma_{Helper::UndefValueDouble};
  double left_sideband_external_{Helper::UndefValueDouble};
  double left_sideband_internal_{Helper::UndefValueDouble};
  double right_sideband_external_{Helper::UndefValueDouble};
  double right_sideband_internal_{Helper::UndefValueDouble};

  std::string peak_shape_{"Gaus"};
  int bg_pol_n_{2};
};
#endif //QA2_SHAPEFITTER_HPP
