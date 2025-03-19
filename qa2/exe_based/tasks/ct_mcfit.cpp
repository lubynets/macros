//
// Created by oleksii on 04.03.25.
//
#include "Helper.hpp"

#include <TCanvas.h>
#include <TLegend.h>
#include <TROOT.h>

#include <iostream>

using namespace Helper;

void ct_mcfit(const std::string& fileNameEff, const std::string& fileNameYield, bool isSaveToRoot) {
  TString currentMacroPath = __FILE__;
  TString directory = currentMacroPath(0, currentMacroPath.Last('/'));
  gROOT->Macro( directory + "/../styles/mc_qa2.dpg.style.cc" );

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
  histoRecCorrected->Fit(fitRecCorrected, "", "", 0.2, 2);

  TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.82);
  leg->AddEntry(histoGen, "MC", "P");
  leg->AddEntry(histoRecCorrected, "Rec (corrected)", "L");
  leg->AddEntry(fitRecCorrected, "Fit to Rec", "L");

  TCanvas cc("cc", "cc", 1200, 800);
  cc.SetLogy();

  histoGen->SetLineColor(kBlack);
  histoGen->SetMarkerStyle(kOpenCircle);
  histoGen->SetMarkerColor(kBlack);
  histoRecCorrected->SetLineColor(kBlue);

  histoGen->GetXaxis()->SetTitle("T (ps)"); // TODO ad. hoc. to be changed in the main code for histo production
  histoGen->Draw("P");
  histoRecCorrected->Draw("same");
  leg->Draw("same");
  AddOneLineText(promptness, {0.69, 0.82, 0.82, 0.90});
  const std::string lifetimeFitRecCorrected = "#tau_{#Lambda_{c}} [Fit] = (" +
                                              to_string_with_significant_figures(fitRecCorrected->GetParameter(1)*1000, 5) +
                                              " #pm " +
                                              to_string_with_significant_figures(fitRecCorrected->GetParError(1)*1000, 2) +
                                              ") fs";
  AddOneLineText(lifetimeFitRecCorrected, {0.69, 0.62, 0.82, 0.70});

  const std::string lifetimeFitMC = "#tau_{#Lambda_{c}} [MC] = (" +
                                    to_string_with_significant_figures(fitMC->GetParameter(1)*1000, 5) +
                                    " #pm " +
                                    to_string_with_significant_figures(fitMC->GetParError(1)*1000, 2) +
                                    ") fs";
  AddOneLineText(lifetimeFitMC, {0.69, 0.57, 0.82, 0.65});

  const std::string lifetimePdg = "#tau_{#Lambda_{c}} [PDG] = (202.6 #pm 1.0) fs";
  AddOneLineText(lifetimePdg, {0.69, 0.52, 0.82, 0.60});

  cc.Print((fileOutName + ".pdf").c_str(), "pdf");

  //=========== build ratio ==============================
  TH1D* histoRecCorrected2RecFit = dynamic_cast<TH1D*>(histoRecCorrected->Clone());
  histoRecCorrected2RecFit->Divide(fitRecCorrected);
  histoRecCorrected2RecFit->GetYaxis()->SetTitle("Ratio");

  TCanvas ccRatio2RecFit("ccRatio2RecFit", "", 1200, 800);
  TF1* oneline = new TF1("oneline", "[0]", 0, 2);
  oneline->SetParameter(0, 1);
  oneline->SetLineColor(kBlack);
  oneline->SetLineStyle(7);
  histoRecCorrected2RecFit->Draw("HIST PE");
  oneline->Draw("same");
  ccRatio2RecFit.Print((fileOutName + "Ratio2RecFit.pdf").c_str(), "pdf");
 // ---------------------------------------------------------
  TH1D* histoMC2MCFit = dynamic_cast<TH1D*>(histoGen->Clone());
  histoMC2MCFit->Divide(fitMC);
  histoMC2MCFit->GetYaxis()->SetTitle("Ratio");
  histoMC2MCFit->GetYaxis()->SetRangeUser(0, 2);

  TH1D* histoRecCorrected2MCFit = dynamic_cast<TH1D*>(histoRecCorrected->Clone());
  histoRecCorrected2MCFit->Divide(fitMC);

  TF1* fitRecCorrected2MCFit = new TF1("fitRecCorrected2MCFit", "[0]*TMath::Exp(-x/[1])", 0.2, 2);
  fitRecCorrected2MCFit->SetParameter(0, fitRecCorrected->GetParameter(0) / fitMC->GetParameter(0));
  fitRecCorrected2MCFit->SetParameter(1, fitMC->GetParameter(1) * fitRecCorrected->GetParameter(1) / (fitMC->GetParameter(1) - fitRecCorrected->GetParameter(1)));
  std::cout << fitRecCorrected2MCFit->GetParameter(0) << "\t" << fitRecCorrected2MCFit->GetParameter(1) << "\n";

  TCanvas ccRatio2MCFit("ccRatio2MCFit", "", 1200, 800);
  histoMC2MCFit->Draw("HIST PE");
  histoRecCorrected2MCFit->Draw("same HIST PE");
  fitRecCorrected2MCFit->Draw("same");
  oneline->Draw("same");

  TLegend* legratio = new TLegend(0.3, 0.20, 0.5, 0.4);
  legratio->AddEntry(histoMC2MCFit, "MC to MC fit", "P");
  legratio->AddEntry(histoRecCorrected2MCFit, "Rec to MC fit", "L");
  legratio->AddEntry(fitRecCorrected2MCFit, "Rec fit MC fit", "L");
  legratio->Draw("same");

  ccRatio2MCFit.Print((fileOutName + "Ratio2MCFit.pdf").c_str(), "pdf");
  //=========== build ratio ==============================
}

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./ct_mcfit fileNameEff fileNameYield (isSaveRoot=false)" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileNameEff = argv[1];
  const std::string fileNameYield = argv[2];
  const bool isSaveToRoot = argc > 3 ? string_to_bool(argv[3]) : false;

  ct_mcfit(fileNameEff, fileNameYield, isSaveToRoot);

  return 0;
}
