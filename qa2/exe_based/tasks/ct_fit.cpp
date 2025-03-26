//
// Created by oleksii on 19.03.25.
//
#include "Helper.hpp"

#include <TCanvas.h>
#include <TROOT.h>

#include <iostream>

using namespace Helper;

void ct_mcfit(const std::string& fileNameYield, const std::string& fileNameEff, const std::string& mcOrData) {
  TString currentMacroPath = __FILE__;
  TString directory = currentMacroPath(0, currentMacroPath.Last('/'));
  gROOT->Macro( directory + "/../styles/mc_qa2.style.cc" );

  TFile* fileInYield = OpenFileWithNullptrCheck(fileNameYield);
  TFile* fileInEff = OpenFileWithNullptrCheck(fileNameEff);

  bool useIntegral{true};
//  bool useIntegral{false};

  const std::string fitIntegralOption = useIntegral ? "I" : "";

  const std::string fileOutNamePrefix = "ct_fit";

  TH1D* histoEff = GetObjectWithNullptrCheck<TH1D>(fileInEff, "hEffPromptT");

  // original yield histogram
  TH1D* histoYield = GetObjectWithNullptrCheck<TH1D>(fileInYield, "hRawYields");
  histoYield->UseCurrentStyle();

  // yield histogram corrected for efficiency
  TH1D* histoYieldCorr4Eff = dynamic_cast<TH1D*>(histoYield->Clone());
  histoYieldCorr4Eff->SetName((static_cast<std::string>(histoYield->GetName()) + "_corr4eff").c_str());
  histoYieldCorr4Eff->GetYaxis()->SetTitle((static_cast<std::string>(histoYieldCorr4Eff->GetYaxis()->GetTitle()) + " corr. eff.").c_str());
  CheckHistogramsForXaxisIdentity(histoYieldCorr4Eff, histoEff);
  for(int iBin=1, nBins = histoYieldCorr4Eff->GetNbinsX(); iBin<=nBins; iBin++) {
    histoYieldCorr4Eff->SetBinContent(iBin, histoYieldCorr4Eff->GetBinContent(iBin) / histoEff->GetBinContent(iBin));
    histoYieldCorr4Eff->SetBinError(iBin, histoYieldCorr4Eff->GetBinError(iBin) / histoEff->GetBinContent(iBin));
  }
  histoYieldCorr4Eff->Scale(100);

  // differential (per x-axis unit) yield histogram
  TH1D* histoYieldDiff = dynamic_cast<TH1D*>(histoYieldCorr4Eff->Clone());
  histoYieldDiff->SetName((static_cast<std::string>(histoYield->GetName()) + "_diff").c_str());
  histoYieldDiff->GetYaxis()->SetTitle("dN/dT (ps^{-1})");
  for(int iBin=1, nBins = histoYieldDiff->GetNbinsX(); iBin<=nBins; iBin++) {
    histoYieldDiff->SetBinContent(iBin, histoYieldDiff->GetBinContent(iBin) / histoYieldDiff->GetBinWidth(iBin));
    histoYieldDiff->SetBinError(iBin, histoYieldDiff->GetBinError(iBin) / histoYieldDiff->GetBinWidth(iBin));
  }

  const double lo = histoYieldDiff->GetBinLowEdge(1) + 1e-3;
  const double hi = histoYieldDiff->GetBinLowEdge(histoYieldDiff->GetNbinsX()+1) - 1e-3;
  auto parEst = EstimateExpoParameters(histoYieldDiff, lo, hi);
  TF1* fitFunc = new TF1("fitFunc", "[0]*TMath::Exp(-x/[1])", lo, hi);
  fitFunc->SetParameters(parEst.first, parEst.second);
  histoYieldDiff->Fit(fitFunc, ("0"+fitIntegralOption).c_str(), "", lo, hi);

  const std::string lifetimeFitValue = "#tau_{#Lambda_{c}} [Fit] = (" +
                                              to_string_with_significant_figures(fitFunc->GetParameter(1)*1000, 3) +
                                              " #pm " +
                                              to_string_with_significant_figures(fitFunc->GetParError(1)*1000, 2) +
                                              ") fs";
  const std::string chi2Value = "#chi^{2} / ndf = " +
                                to_string_with_significant_figures(fitFunc->GetChisquare(), 3) +
                                " / " +
                                std::to_string(fitFunc->GetNDF());

  const std::string lifetimePdg = "#tau_{#Lambda_{c}} [PDG] = (202.6 #pm 1.0) fs";

  TCanvas ccYield("ccYield", "");
  ccYield.SetCanvasSize(1200, 800);
  histoYield->Draw();
  ccYield.Print((fileOutNamePrefix + ".pdf(").c_str(), "pdf");

  TCanvas ccYieldCorr4Eff("ccYieldCorr4Eff", "");
  ccYieldCorr4Eff.SetCanvasSize(1200, 800);
  histoYieldCorr4Eff->Draw();
  ccYieldCorr4Eff.Print((fileOutNamePrefix + ".pdf").c_str(), "pdf");

  TCanvas ccYieldDiff("ccYieldDiff", "");
  ccYieldDiff.SetLogy();
  ccYieldDiff.SetCanvasSize(1200, 800);
  histoYieldDiff->Draw();
  fitFunc->Draw("same");
  AddMultiLineText({mcOrData, chi2Value, lifetimeFitValue, lifetimePdg}, {0.69, 0.75, 0.82, 0.80});
  ccYieldDiff.Print((fileOutNamePrefix + ".pdf").c_str(), "pdf");

  // ============== build ratio ==============================
  TH1D* histoYieldDiff2Fit = dynamic_cast<TH1D*>(histoYieldDiff->Clone());
  ScalePlotVertically(histoYieldDiff2Fit, histoYieldDiff, 2);
  if(!useIntegral) {
    histoYieldDiff2Fit->Divide(fitFunc);
  } else {
    for(int iBin=1, nBins=histoYieldDiff2Fit->GetNbinsX(); iBin<=nBins; iBin++) {
      const double histoValue = histoYieldDiff2Fit->GetBinContent(iBin);
      const double histoError = histoYieldDiff2Fit->GetBinError(iBin);
      const double lo = histoYieldDiff2Fit->GetBinLowEdge(iBin);
      const double hi = histoYieldDiff2Fit->GetBinLowEdge(iBin+1);
      const double funcAverage = fitFunc->Integral(lo, hi) / (hi-lo);
      histoYieldDiff2Fit->SetBinContent(iBin, histoValue/funcAverage);
      histoYieldDiff2Fit->SetBinError(iBin, histoError/funcAverage);
    }
  }
  histoYieldDiff2Fit->GetYaxis()->SetTitle("Data / Fit");

  TCanvas ccRatio("ccRatio", "");
  ScaleCanvasVertically(&ccRatio, &ccYieldDiff, 2);
  TF1 oneline("oneline", "[0]", 0, 2);
  oneline.SetParameter(0, 1);
  oneline.SetLineColor(kBlack);
  oneline.SetLineStyle(7);
  histoYieldDiff2Fit->Draw("HIST PE");
  oneline.Draw("same");
  ccRatio.Print((fileOutNamePrefix + ".pdf)").c_str(), "pdf");
  // ============== end build ratio ==========================

  delete fitFunc;
}

int main(int argc, char* argv[]) {
  if (argc < 4) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./ct_mcfit fileNameYield fileNameEff mcOrData" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileNameYield = argv[1];
  const std::string fileNameEff = argv[2];
  const std::string mcOrData = argv[3];

  ct_mcfit(fileNameYield, fileNameEff, mcOrData);

  return 0;
}