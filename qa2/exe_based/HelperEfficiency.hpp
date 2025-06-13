//
// Created by oleksii on 13.06.25.
//

#ifndef QA2_HELPEREFFICIENCY_HPP
#define QA2_HELPEREFFICIENCY_HPP

#include "Helper.hpp"

#include <TH1.h>

namespace HelperEfficiency {

inline std::pair<TH1*, TH1*> EvaluateEfficiencyHisto(TH1* hNum, TH1* hDen) {
  const int nBins = hNum->GetNbinsX();
  Helper::CheckHistogramsForXaxisIdentity(hNum, hDen);

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

  hEff->GetYaxis()->SetTitle("#varepsilon, %");
  hEff->Scale(100);
  hEff->SetTitle("");

  hRelErr->GetYaxis()->SetTitle("#varepsilon_{#varepsilon}");
  hRelErr->SetTitle("");

  return std::make_pair(hEff, hRelErr);
}

}

#endif //QA2_HELPEREFFICIENCY_HPP
