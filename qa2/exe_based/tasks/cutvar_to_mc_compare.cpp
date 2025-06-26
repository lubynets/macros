//
// Created by oleksii on 26.06.25.
//
#include "Helper.hpp"

#include <TFile.h>
#include <TLegend.h>

#include <iostream>
#include <string>

using namespace Helper;

void cutvar_to_mc_compare(const std::string& fileNameCutVar, const std::string& fileNameMC) {
  LoadMacro("styles/mc_qa2.style.cc");

  const std::vector<double> lifetimeRanges = {0.2, 0.35, 0.5, 0.7, 0.9, 1.6};

  TFile* fileCutVar = OpenFileWithNullptrCheck(fileNameCutVar);
  TFile* fileMC = OpenFileWithNullptrCheck(fileNameMC);

  TH1* hPromptCutVar = GetObjectWithNullptrCheck<TH1>(fileCutVar, "hCorrYieldsPrompt");
  TH1* hNonPromptCutVar = GetObjectWithNullptrCheck<TH1>(fileCutVar, "hCorrYieldsNonPrompt");
  for(auto& h : {hPromptCutVar, hNonPromptCutVar}) {
    h->UseCurrentStyle();
    h->SetLineColor(kRed);
    h->SetMarkerColor(kRed);
    h->GetYaxis()->SetMoreLogLabels();
  }

  TH1* hPromptMC = GetObjectWithNullptrCheck<TH1>(fileMC, "gen/prompt/hT");
  TH1* hNonPromptMC = GetObjectWithNullptrCheck<TH1>(fileMC, "gen/nonprompt/hT");
  for(auto& h : {hPromptMC, hNonPromptMC}) {
    h->UseCurrentStyle();
    h->SetLineColor(kBlue);
    h->SetMarkerColor(kBlue);
    h->GetYaxis()->SetMoreLogLabels();
  }
  hPromptMC = dynamic_cast<TH1D*>(hPromptMC->Rebin(lifetimeRanges.size() - 1,hPromptMC->GetName(),lifetimeRanges.data()));
  hNonPromptMC = dynamic_cast<TH1D*>(hNonPromptMC->Rebin(lifetimeRanges.size() - 1,hNonPromptMC->GetName(),lifetimeRanges.data()));

  CustomizeHistogramsYRange({hPromptCutVar, hPromptMC}, true);
  CustomizeHistogramsYRange({hNonPromptCutVar, hNonPromptMC}, true);

  TCanvas ccPrompt("ccPrompt", "");
  ccPrompt.SetLogy();
  ccPrompt.SetCanvasSize(1200, 800);
  hPromptMC->Draw("E1");
  hPromptCutVar->Draw("E1 same");
  TLegend legPrompt(0.7, 0.7, 0.9, 0.9);
  legPrompt.SetHeader("prompt");
  legPrompt.AddEntry(hPromptMC, "MC", "PL");
  legPrompt.AddEntry(hPromptCutVar, "Cut Var", "PL");
  legPrompt.Draw("same");
  ccPrompt.Print("compar.pdf(", "pdf");

  TCanvas ccNonPrompt("ccNonPrompt", "");
  ccNonPrompt.SetLogy();
  ccNonPrompt.SetCanvasSize(1200, 800);
  hNonPromptMC->Draw("E1");
  hNonPromptCutVar->Draw("E1 same");
  TLegend legNonPrompt(0.7, 0.7, 0.9, 0.9);
  legNonPrompt.SetHeader("nonprompt");
  legNonPrompt.AddEntry(hNonPromptMC, "MC", "PL");
  legNonPrompt.AddEntry(hNonPromptCutVar, "Cut Var", "PL");
  legNonPrompt.Draw("same");
  ccNonPrompt.Print("compar.pdf)", "pdf");
}

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./cutvar_to_mc_compare fileNameCutVar fileNameMC" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileNameCutVar = argv[1];
  const std::string fileNameMC = argv[2];

  cutvar_to_mc_compare(fileNameCutVar, fileNameMC);

  return 0;
}