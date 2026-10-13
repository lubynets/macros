//
// Created by oleksii on 13.10.2026
//

#include "ForwardFoldedExpo.hpp"

#include <cmath>

void ForwardFoldedExpo::Init() {
  n_bins_ = response_matrix_->GetNbinsX();
}

double ForwardFoldedExpo::IntegrateExpo(double from, double to) {
  return expo_par_A_ * expo_par_tau_ * (std::exp(-from/expo_par_tau_) - std::exp(-to/expo_par_tau_));
}

double ForwardFoldedExpo::IntegrateExpo(int iBin) {
  return IntegrateExpo(response_matrix_->GetYaxis()->GetBinLowEdge(iBin), response_matrix_->GetYaxis()->GetBinLowEdge(iBin+1));
}

double ForwardFoldedExpo::operator()(double* x, double* par) {
  expo_par_A_ = par[0];
  expo_par_tau_ = par[1];
  const int binContaningX = response_matrix_->GetXaxis()->FindBin(*x);
  double result{};
  for(int iBin=1; iBin<=n_bins_; ++iBin) {
    result += response_matrix_->GetBinContent(binContaningX, iBin) * IntegrateExpo(iBin);
  }

  return result;
}
