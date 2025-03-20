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

  const std::vector<std::string> weightUsages {
    "wo_weight",
    "with_weight"
  };

  // ad hoc
  if(!((fileNameEff.find(".Full.") != std::string::npos && mcOrData == "Data") || (fileNameEff.find(".I.") != std::string::npos && mcOrData == "Mc"))) {
    std::cout << "WARNING! Something wrong with input efficiency histogram and type of input (mc or data)\n";
  }

  const std::string fileOutNamePrefix = "ct_fit_";

  TH1D* histoEff = GetObjectWithNullptrCheck<TH1D>(fileInEff, "hEffPromptT");

  for(auto& wu : weightUsages) {
    const std::string weightsApplication = wu == "wo_weight" ? "coarse weights" : "fine weights";

    // original yield histogram
    TH1D* histoYield = GetObjectWithNullptrCheck<TH1D>(fileInYield, ("histoYieldSignal_" + wu).c_str());

    // yield histogram corrected for efficiency (if weights not used, otherwise it is already corrected for efficiency)
    TH1D* histoYieldCorr4Eff = dynamic_cast<TH1D*>(histoYield->Clone());
    histoYieldCorr4Eff->SetName((static_cast<std::string>(histoYield->GetName()) + "_corr4eff").c_str());
    if(wu == "wo_weight") {
      CheckHistogramsForXaxisIdentity(histoYieldCorr4Eff, histoEff);
      for(int iBin=1, nBins = histoYieldCorr4Eff->GetNbinsX(); iBin<=nBins; iBin++) {
        histoYieldCorr4Eff->SetBinContent(iBin, histoYieldCorr4Eff->GetBinContent(iBin) / histoEff->GetBinContent(iBin));
        histoYieldCorr4Eff->SetBinError(iBin, histoYieldCorr4Eff->GetBinError(iBin) / histoEff->GetBinContent(iBin));
      }
      histoYieldCorr4Eff->Scale(100);
    }

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
    histoYieldDiff->Fit(fitFunc, "0I", "", lo, hi);

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
    ccYield.Print((fileOutNamePrefix + wu + ".pdf(").c_str(), "pdf");

    TCanvas ccYieldCorr4Eff("ccYieldCorr4Eff", "");
    ccYieldCorr4Eff.SetCanvasSize(1200, 800);
    histoYieldCorr4Eff->Draw();
    ccYieldCorr4Eff.Print((fileOutNamePrefix + wu + ".pdf").c_str(), "pdf");

    TCanvas ccYieldDiff("ccYieldDiff", "");
    ccYieldDiff.SetLogy();
    ccYieldDiff.SetCanvasSize(1200, 800);
    histoYieldDiff->Draw();
    fitFunc->Draw("same");
    const float oltX1 = 0.69, oltX2 = 0.82, oltY2 = 0.80, oltYstep = 0.05;
    AddOneLineText(mcOrData, {oltX1, oltY2-0*oltYstep, oltX2, oltY2-1*oltYstep});
    AddOneLineText(weightsApplication, {oltX1, oltY2-1*oltYstep, oltX2, oltY2-2*oltYstep});
    AddOneLineText(chi2Value, {oltX1, oltY2-2*oltYstep, oltX2, oltY2-3*oltYstep});
    AddOneLineText(lifetimeFitValue, {oltX1, oltY2-3*oltYstep, oltX2, oltY2-4*oltYstep});
    AddOneLineText(lifetimePdg, {oltX1, oltY2-4*oltYstep, oltX2, oltY2-5*oltYstep});
    ccYieldDiff.Print((fileOutNamePrefix + wu + ".pdf)").c_str(), "pdf");
  }
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