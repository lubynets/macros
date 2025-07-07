//
// Created by oleksii on 26.06.25.
//
#include "Helper.hpp"

#include <TFile.h>
#include <TLegend.h>

#include <iostream>
#include <string>

using namespace Helper;

void corrected_yields_qa2(const std::string& fileNameCutVar, const std::string& fileNameMC) {
  LoadMacro("styles/mc_qa2.style.cc");

  const bool isMc = !fileNameMC.empty();

  const std::vector<double> lifetimeRanges = {0.2, 0.35, 0.5, 0.7, 0.9, 1.6};

  TFile* fileCutVar = OpenFileWithNullptrCheck(fileNameCutVar);
  TFile* fileMC = isMc ? OpenFileWithNullptrCheck(fileNameMC) : nullptr;

  TH1* hPromptCutVar = GetObjectWithNullptrCheck<TH1>(fileCutVar, "hCorrYieldsPrompt");
  TH1* hNonPromptCutVar = GetObjectWithNullptrCheck<TH1>(fileCutVar, "hCorrYieldsNonPrompt");
  for(auto& h : {hPromptCutVar, hNonPromptCutVar}) {
    h->UseCurrentStyle();
    h->SetLineColor(kRed);
    h->SetMarkerColor(kRed);
    h->GetYaxis()->SetMoreLogLabels();
  }

  TH1* hPromptMC = isMc ? GetObjectWithNullptrCheck<TH1>(fileMC, "gen/prompt/hT") : nullptr;
  TH1* hNonPromptMC = isMc ? GetObjectWithNullptrCheck<TH1>(fileMC, "gen/nonprompt/hT") : nullptr;
  if(isMc) {
    for (auto& h: {hPromptMC, hNonPromptMC}) {
      h->UseCurrentStyle();
      h->SetLineColor(kBlue);
      h->SetMarkerColor(kBlue);
      h->GetYaxis()->SetMoreLogLabels();
    } // {hPromptMC, hNonPromptMC}
    hPromptMC = dynamic_cast<TH1D*>(hPromptMC->Rebin(lifetimeRanges.size() - 1,hPromptMC->GetName(),lifetimeRanges.data()));
    hNonPromptMC = dynamic_cast<TH1D*>(hNonPromptMC->Rebin(lifetimeRanges.size() - 1,hNonPromptMC->GetName(),lifetimeRanges.data()));

    CustomizeHistogramsYRange({hPromptCutVar, hPromptMC}, true);
    CustomizeHistogramsYRange({hNonPromptCutVar, hNonPromptMC}, true);
  } // isMc

  auto PrintCorrectedYields = [&](TH1* hCutVar, TH1* hMc, const std::string& promptness, const std::string& priBra) {
    TCanvas cc("cc", "");
    cc.SetLogy();
    cc.SetCanvasSize(1200, 800);
    if(isMc) hMc->Draw("E1");
    hCutVar->Draw("E1 same");
    TLegend legPrompt(0.7, 0.7, 0.9, 0.9);
    legPrompt.SetHeader(promptness.c_str());
    if(isMc) legPrompt.AddEntry(hMc, "MC", "PL");
    legPrompt.AddEntry(hCutVar, "Cut Var", "PL");
    legPrompt.Draw("same");
    cc.Print(("compar.pdf" + priBra).c_str(), "pdf");
  };

  PrintCorrectedYields(hPromptCutVar, hPromptMC, "prompt", "(");
  PrintCorrectedYields(hNonPromptCutVar, hNonPromptMC, "nonprompt", ")");
}

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./corrected_yields_qa2 fileNameCutVar (fileNameMC="")" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileNameCutVar = argv[1];
  const std::string fileNameMC = argv[2];

  corrected_yields_qa2(fileNameCutVar, fileNameMC);

  return 0;
}