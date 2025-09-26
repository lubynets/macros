//
// Created by oleksii on 03.07.25.
//
#include "HelperGeneral.hpp"
#include "HelperMath.hpp"
#include "HelperPlot.hpp"

#include <TGraphErrors.h>
#include <TH1.h>
#include <TStyle.h>

#include <iostream>
#include <numeric>

using namespace HelperGeneral;
using namespace HelperMath;
using namespace HelperPlot;

double EvaluateAverageExcludingOutliers(const std::vector<double>& values, const std::vector<double>& errors, double chi2Max = 1.);
double EvaluateAverageExcludingOutliers(const TGraphErrors* graph, double from=-1e9, double to=1e9, double chi2Max=1.);

void complex_vs_bdt_pdfer(const std::string& fileNameTemplate, const std::string& targetSignal, const std::string& wise){
  LoadMacro("styles/mc_qa2.style.cc");
  gStyle->SetMarkerSize(0.6);
  gStyle->SetNdivisions(315, "X");
  gStyle->SetLineWidth(1);
  gStyle->SetPadGridX(true);
  //=================================================================
  std::vector<double> bdtScores;
  for(int i=0; i<=99; i++) {
    bdtScores.emplace_back(0.01*i);
  }
  //=================================================================

  const std::vector<std::string> variables {
    "RawYieldsSignal",
//    "RawYieldsSignalCounted",
//    "RawYieldsSigma",
//    "RawYieldsMean",
//    "RawYieldsSignificance",
//    "RawYieldsSgnOverBkg",
//    "RawYieldsBkg",
//    "RawYieldsChiSquareBkg",
//    "RawYieldsChiSquareTotal"
  };

  TFile* fileMarkup = OpenFileWithNullptrCheck(fileNameTemplate + "." + targetSignal + "gt" + to_string_with_precision(bdtScores.at(0), 2) + ".root");
  TH1* histoMarkup = GetObjectWithNullptrCheck<TH1>(fileMarkup, "h" + variables.at(0));
  std::vector<float> lifeTimeRanges;
  for(int iBin=1; iBin<=histoMarkup->GetNbinsX()+1; ++iBin) {
    lifeTimeRanges.emplace_back(histoMarkup->GetBinLowEdge(iBin));
  }

  HelperMath::tensor2<TGraphErrors*> graphVar = HelperMath::make_tensor<TGraphErrors*, 2>({variables.size(), lifeTimeRanges.size()-1}, nullptr);
  for(int iVar=0, nVars=variables.size(); iVar<nVars; ++iVar) {
    histoMarkup = GetObjectWithNullptrCheck<TH1>(fileMarkup, "h" + variables.at(iVar));
    const std::string yAxisTitle = histoMarkup->GetYaxis()->GetTitle();
    for(int iT=0, nTs=lifeTimeRanges.size()-1; iT<nTs; ++iT) {
      graphVar.at(iVar).at(iT) = new TGraphErrors();
      auto gra = graphVar.at(iVar).at(iT);
      gra->SetName(("gr" + variables.at(iVar) + "_" + + "_T" + std::to_string(iT)).c_str());
      gra->SetTitle(("bin #" + std::to_string(iT+1) + "#; T#in (" + to_string_with_precision(lifeTimeRanges.at(iT), 2) + "#; " + to_string_with_precision(lifeTimeRanges.at(iT+1), 2) + ") ps").c_str());
      gra->GetXaxis()->SetTitle(("bdt score " + targetSignal).c_str());
      gra->GetYaxis()->SetTitle(yAxisTitle.c_str());
    } // lifeTimeRanges
  } // variables
  fileMarkup->Close();

  for(const auto& score : bdtScores) {
    TFile* fileIn = OpenFileWithNullptrCheck(fileNameTemplate + "." + targetSignal + "gt" + to_string_with_precision(score, 2) + ".root");
    for(int iVar=0, nVars=variables.size(); iVar<nVars; ++iVar) {
      TH1* histoVar = GetObjectWithNullptrCheck<TH1>(fileIn, "h" + variables.at(iVar));
      for(int iT=0, nTs=lifeTimeRanges.size()-1; iT<nTs; ++iT) {
        auto gra = graphVar.at(iVar).at(iT);
        gra->SetPoint(gra->GetN(), score, histoVar->GetBinContent(iT+1));
        gra->SetPointError(gra->GetN()-1, 0, histoVar->GetBinError(iT+1));
      } // lifetimeRanges
    } // variables
    fileIn->Close();
  } // bdtScores

  std::string priBra, ccName;
  TH1F* hMeanFix = new TH1F("hRawYieldsMean", "hRawYieldsMean", lifeTimeRanges.size()-1, lifeTimeRanges.data());
  TH1F* hSigmaFix = new TH1F("hRawYieldsSigma", "hRawYieldsSigma", lifeTimeRanges.size()-1, lifeTimeRanges.data());
  for(int iVar=0, nVars=variables.size(); iVar<nVars; ++iVar) {
    for(int iT=0, nTs=lifeTimeRanges.size()-1; iT<nTs; ++iT) {
      if(wise == "var") {
        priBra = lifeTimeRanges.size()-1 == 1 ? "" : iT == 0 ? "(" : iT == lifeTimeRanges.size()-2 ? ")" : "";
        ccName = variables.at(iVar) + "_vs_" + targetSignal;
      } else {
        priBra = variables.size() == 1 ? "" : iVar == 0 ? "(" : iVar == variables.size()-1 ? ")" : "";
        ccName = "T_" + to_string_with_precision(lifeTimeRanges.at(iT), 2) + "_vs_" + targetSignal;
      }
      TCanvas cc("cc", "");
      cc.SetCanvasSize(1200, 800);
      const auto gra = graphVar.at(iVar).at(iT);
      gra->Draw("APE");
      if(variables.at(iVar) == "RawYieldsMean") HorizontalLine4Graph(massLambdaC, gra)->Draw("same");
      if(variables.at(iVar) == "RawYieldsMean" || variables.at(iVar) == "RawYieldsSigma") {
        const double ave = EvaluateAverageExcludingOutliers(gra, -0.01, 0.101);
        auto aveLine = HorizontalLine4Graph(ave, gra);
        aveLine->SetLineColor(kBlack);
        aveLine->SetLineStyle(7);
        aveLine->Draw("same");
        auto hFix = variables.at(iVar) == "RawYieldsMean" ? hMeanFix : hSigmaFix;
        hFix->SetBinContent(iT+1, ave);
      }
      cc.Print((ccName + ".pdf" + priBra).c_str(), "pdf");
    } // lifeTimeRanges
  } // variables
  TFile* fileOutFix = TFile::Open((fileNameTemplate + "." + targetSignal + "_fix.root").c_str(), "recreate");
  hMeanFix->Write();
  hSigmaFix->Write();
  fileOutFix->Close();
  delete hMeanFix;
  delete hSigmaFix;
}

double EvaluateAverageExcludingOutliers(const std::vector<double>& values, const std::vector<double>& errors, double chi2Max) {
  if(values.size() != errors.size()) throw std::runtime_error("EvaluateAverageExcludingOutliers() - values.size() != errors.size()");
  const double averagePreliminary = std::accumulate(values.begin(), values.end(), 0.) / values.size();
  double sum{0.};
  int count{0};
  for(int iEl=0, nEls=values.size(); iEl<nEls; ++iEl) {
    if(std::fabs(values.at(iEl) - averagePreliminary) / errors.at(iEl) > std::sqrt(chi2Max)) continue;
    sum += values.at(iEl);
    ++count;
  }

  return count > 0 ? sum/count : averagePreliminary;
}

double EvaluateAverageExcludingOutliers(const TGraphErrors* graph, double from, double to, double chi2Max) {
  std::vector<double> values, errors;
  for(int iPoint=0, nPoints=graph->GetN(); iPoint<nPoints; ++iPoint) {
    const double x = graph->GetPointX(iPoint);
    if(x<from || x>to) continue;
    values.emplace_back(graph->GetPointY(iPoint));
    errors.emplace_back(graph->GetErrorY(iPoint));
  }

  return EvaluateAverageExcludingOutliers(values, errors, chi2Max);
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./complex_vs_bdt_pdfer fileNameTemplate (targetSignal=NP wise=var)" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileNameTemplate = argv[1];
  const std::string targetSignal = argc > 2 ? argv[2] : "NP";
  const std::string wise = argc > 3 ? argv[3] : "var";
  if(targetSignal != "P" && targetSignal != "NP") throw std::runtime_error("main(): targetSignal must be either 'P' or 'NP'");
  if(wise != "var" && wise != "T") throw std::runtime_error("main(): wise must be either 'var' or 'T'");

  complex_vs_bdt_pdfer(fileNameTemplate, targetSignal, wise);

  return 0;
}