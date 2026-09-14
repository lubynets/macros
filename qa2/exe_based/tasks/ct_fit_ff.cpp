//
// Created by oleksii on 14.09.2026
//
#include "HelperGeneral.hpp"

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

#include <cmath>
#include <iostream>
#include <string>

using namespace HelperGeneral;

class ForwardFoldedExpo {
public:
  void SetResponseMatrix(TH2* rm) { response_matrix_ = rm; }
  void Init();

  double operator()(double* x, double* par);

private:
  double IntegrateExpo(double from, double to);
  double IntegrateExpo(int iBin);

  int n_bins_{};
  double expo_par_A_{};
  double expo_par_tau_{};
  TH2* response_matrix_{};
};

void ForwardFoldedExpo::Init() {
  n_bins_ = response_matrix_->GetNbinsX();
}

double ForwardFoldedExpo::IntegrateExpo(double from, double to) {
  return expo_par_A_ * expo_par_tau_ * (std::exp(-from/expo_par_tau_) - std::exp(-to/expo_par_tau_));
}

double ForwardFoldedExpo::IntegrateExpo(int iBin) {
  return IntegrateExpo(response_matrix_->GetXaxis()->GetBinLowEdge(iBin), response_matrix_->GetXaxis()->GetBinLowEdge(iBin+1));
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

void ct_fit_ff(const std::string& fileNameYield, const std::string& fileNameResponseMatrix, const std::string& histoNameResponseMatrix) {
  TFile* fileYield = OpenFileWithNullptrCheck(fileNameYield);
  TH1* histoYield = GetObjectWithNullptrCheck<TH1>(fileYield, "hCorrYieldsPrompt");

  TFile* fileRespMatrix = OpenFileWithNullptrCheck(fileNameResponseMatrix);
  TH2* histoRespMatrix = GetObjectWithNullptrCheck<TH2>(fileRespMatrix, histoNameResponseMatrix);

  CheckHistogramsForAxisIdentity<TH2, TH2>(histoRespMatrix, nullptr, "XY");
  CheckHistogramsForAxisIdentity(histoRespMatrix, histoYield, "X");

//   TH1* histEff = histoRespMatrix->ProjectionX(); // for expo parameters estimate only


  fileRespMatrix->Close();
  fileYield->Close();
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
