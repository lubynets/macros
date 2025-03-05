//
// Created by oleksii on 04.03.25.
//
#include "Helper.hpp"

#include <TCanvas.h>
#include <TLegend.h>
#include <TROOT.h>

#include <iostream>

using namespace Helper;

std::pair<float, float> EstimateExpoParameters(TH1* h, float lo, float hi);

void fit_yield(const std::string& fileNameEff, const std::string& fileNameYield, bool isSaveToRoot) {
  TString currentMacroPath = __FILE__;
  TString directory = currentMacroPath(0, currentMacroPath.Last('/'));
  gROOT->Macro( directory + "/../styles/mc_qa2.style.cc" );

  TFile* fileInEff = TFile::Open(fileNameEff.c_str());
  TFile* fileInYield = TFile::Open(fileNameYield.c_str());
  if(fileInEff == nullptr || fileInYield == nullptr) {
    throw std::runtime_error("fileInEff == nullptr || fileInYield == nullptr");
  }

  const std::string promptness = "Prompt";

  const std::string fileOutName = "yield_fit." + promptness;

  TH1D* histoRec = fileInYield->Get<TH1D>(("rec_T_" + promptness).c_str());
  TH1D* histoGen = fileInYield->Get<TH1D>(("gen_T_" + promptness).c_str());
  TH1D* histoEff = fileInEff->Get<TH1D>(("hEff" + promptness + "T").c_str());

  for(auto& histo : {histoRec, histoGen, histoEff}) {
    histo->Sumw2();
    histo->SetLineWidth(3);
  }

  TH1D* histoRecCorrected = (TH1D*)histoRec->Clone();
  histoRecCorrected->Divide(histoEff);
  histoRecCorrected->Scale(100);

  TF1* fitMC = new TF1("fitMC", "[0]*TMath::Exp(-x/[1])", 0, 2);
  auto pars_est_MC = EstimateExpoParameters(histoRecCorrected, 0.2, 2);
  fitMC->SetParameters(pars_est_MC.first, pars_est_MC.second);
  histoGen->Fit(fitMC, "0", "", 0.2, 2);

  TF1* fitRecCorrected = new TF1("fitRecCorrected", "[0]*TMath::Exp(-x/[1])", 0, 2);
  auto pars_est_Rec = EstimateExpoParameters(histoRecCorrected, 0.2, 2);
  fitRecCorrected->SetParameters(pars_est_Rec.first, pars_est_Rec.second);
  histoRecCorrected->Fit(fitRecCorrected, "0", "", 0.2, 2);

  TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.82);
  leg->AddEntry(histoGen, "MC", "L");
  leg->AddEntry(histoRecCorrected, "Rec (corrected)", "L");
  leg->AddEntry(fitRecCorrected, "Fit to Rec", "L");

  TCanvas cc("cc", "cc", 1200, 800);
  cc.SetLogy();

  histoGen->SetLineColor(kBlack);
  histoRecCorrected->SetLineColor(kBlue);

  histoGen->Draw("");
  histoRecCorrected->Draw("same");
  fitRecCorrected->Draw("same");
  leg->Draw("same");
  AddOneLineText(promptness, 0.69, 0.82, 0.82, 0.90);
  const std::string lifetimeFitRecCorrected = "#tau_{#Lambda_{c}} [Fit] = (" +
                                              to_string_with_significant_figures(fitRecCorrected->GetParameter(1)*1000, 5) +
                                              " #pm " +
                                              to_string_with_significant_figures(fitRecCorrected->GetParError(1)*1000, 2) +
                                              ") fs";
  AddOneLineText(lifetimeFitRecCorrected, 0.69, 0.62, 0.82, 0.70);

  const std::string lifetimeFitMC = "#tau_{#Lambda_{c}} [MC] = (" +
                                    to_string_with_significant_figures(fitMC->GetParameter(1)*1000, 5) +
                                    " #pm " +
                                    to_string_with_significant_figures(fitMC->GetParError(1)*1000, 2) +
                                    ") fs";
  AddOneLineText(lifetimeFitMC, 0.69, 0.57, 0.82, 0.65);

  const std::string lifetimePdg = "#tau_{#Lambda_{c}} [PDG] = (202.6 #pm 1.0) fs";
  AddOneLineText(lifetimePdg, 0.69, 0.52, 0.82, 0.60);

  cc.Print((fileOutName + ".pdf").c_str(), "pdf");
}

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./fit_yield fileNameEff fileNameYield (isSaveRoot=false)" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileNameEff = argv[1];
  const std::string fileNameYield = argv[2];
  const bool isSaveToRoot = argc > 3 ? string_to_bool(argv[3]) : false;

  fit_yield(fileNameEff, fileNameYield, isSaveToRoot);

  return 0;
}

std::pair<float, float> EstimateExpoParameters(TH1* h, float lo, float hi) {
  const int ilo = h->FindBin(lo);
  const int ihi = h->FindBin(hi);
  const float flo = h->GetBinContent(ilo)/* * h->GetBinWidth(ilo)*/;
  const float fhi = h->GetBinContent(ihi)/* * h->GetBinWidth(ihi)*/;
  const float tau = (hi-lo)/std::log(flo/fhi);
  const float A = flo / std::exp(-lo/tau);
  return std::make_pair(A, tau);
}