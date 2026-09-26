//
// Created by oleksii on 14.09.2026
//
#include "HelperGeneral.hpp"

#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace HelperGeneral;

void ExtendHistoWithEmptyBinsLeft(TH1*& histo, const std::vector<double>& additionalBinEdges);

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

[[deprecated]] void zeroNonDiagonalBins(TH2* histo) {
  const int nBins = histo->GetNbinsX();
  for(int iBin=1; iBin<=nBins; ++iBin) {
    for(int jBin=1; jBin<=nBins; ++jBin) {
      if(iBin != jBin) histo->SetBinContent(iBin, jBin, 0.);
    }
  }
}

void ct_fit_ff(const std::string& fileNameYield, const std::string& fileNameResponseMatrix, const std::string& histoNameResponseMatrix) {
  TFile* fileYield = OpenFileWithNullptrCheck(fileNameYield);
  TH1* histoYield = GetObjectWithNullptrCheck<TH1>(fileYield, "hCorrYieldsPrompt");
//   ExtendHistoWithEmptyBinsLeft(histoYield, {0., 0.2});

  histoYield->SetMarkerColor(kBlue);
  histoYield->SetLineColor(kBlue);
  histoYield->SetMarkerStyle(kFullSquare);
  histoYield->SetMarkerSize(1.6);
  histoYield->GetYaxis()->SetTitle("semi-corrected yield prompt");

  TFile* fileRespMatrix = OpenFileWithNullptrCheck(fileNameResponseMatrix);
  TH2* histoRespMatrix = GetObjectWithNullptrCheck<TH2>(fileRespMatrix, histoNameResponseMatrix);
//   zeroNonDiagonalBins(histoRespMatrix);

  CheckHistogramsForAxisIdentity<TH2, TH2>(histoRespMatrix, nullptr, "XY");
  CheckHistogramsForAxisIdentity(histoRespMatrix, histoYield, "X");

  const double lo = 0.2;
  const double hi = 1.8;

  // ------------- expo parameters estimate --------------------------------------------
  TH1* histEff = histoRespMatrix->ProjectionX(); // for expo parameters estimate only
  TH1* histPreliminaryCorrectedYield = dynamic_cast<TH1*>(histoYield->Clone());
  histPreliminaryCorrectedYield->SetDirectory(nullptr);
  histPreliminaryCorrectedYield->Divide(histEff);
  const double preliminaryYield = histPreliminaryCorrectedYield->Integral();
  const double preliminaryA = preliminaryYield / LifetimeLambdaC / (std::exp(-lo / LifetimeLambdaC) - std::exp(-hi / LifetimeLambdaC));
  // -----------------------------------------------------------------------------------

  ForwardFoldedExpo ffe{};
  ffe.SetResponseMatrix(histoRespMatrix);
  ffe.Init();

  TF1* fitFunc = new TF1("fitFunc", ffe, 0., 2., 2);
  fitFunc->SetParameters(preliminaryA, LifetimeLambdaC);
  fitFunc->SetNpx(1000);

  histoYield->Fit(fitFunc, "", "", lo, hi);

  histoYield->SaveAs("h1.root");

//   //----------debug only-------------------------- TODO remove
//   TFile* fileChi2 = TFile::Open("ct_fit_ff.chi2.root", "update");
//   TH1* hChi2 = GetObjectWithNullptrCheck<TH1>(fileChi2, "hChi2");
//   std::cout << hChi2->GetEntries() << "\t";
//   hChi2->Fill(fitFunc->GetChisquare());
//   std::cout << hChi2->GetEntries() << "\n";
//   fileChi2->WriteObject(hChi2, "hChi2", "Overwrite");
//   fileChi2->Close();
//   //----------------------------------------------------------

  fileRespMatrix->Close();
  fileYield->Close();
}

void ExtendHistoWithEmptyBinsLeft(TH1*& histo, const std::vector<double>& additionalBinEdges) {
  if(additionalBinEdges.size() < 2) {
    throw std::runtime_error("ExtendHistoWithEmptyBinsLeft() - need at least 2 bin edges");
  }
  if(!EqualFloating(additionalBinEdges.back(), histo->GetXaxis()->GetXmin())) {
    throw std::runtime_error("ExtendHistoWithEmptyBinsLeft() - additionalBinEdges.back() != histo->GetXaxis()->GetXmin()");
  }
  if(!std::is_sorted(additionalBinEdges.begin(), additionalBinEdges.end())) {
    throw std::runtime_error("ExtendHistoWithEmptyBinsLeft() - additionalBinEdges is not sorted");
  }
  const int nAddBins = additionalBinEdges.size() - 1;
  std::vector<double> newHistoEdges{additionalBinEdges};
  for(int iBin=1, nBins=histo->GetNbinsX(); iBin<=nBins; ++iBin) {
    newHistoEdges.push_back(histo->GetBinLowEdge(iBin+1));
  }
  TH1* hResult = new TH1D(histo->GetName(), histo->GetTitle(), newHistoEdges.size()-1, newHistoEdges.data());
  hResult->GetXaxis()->SetTitle(histo->GetXaxis()->GetTitle());
  hResult->GetYaxis()->SetTitle(histo->GetYaxis()->GetTitle());
  const bool hasSumw2 = histo->GetSumw2N() > 0;
  if(hasSumw2) hResult->Sumw2();
  for(int iBin=1, nBins=histo->GetNbinsX(); iBin<=nBins; ++iBin) {
    hResult->SetBinContent(iBin+nAddBins, histo->GetBinContent(iBin));
    if(hasSumw2) hResult->SetBinError(iBin+nAddBins, histo->GetBinError(iBin));
  }
  histo = hResult;
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
