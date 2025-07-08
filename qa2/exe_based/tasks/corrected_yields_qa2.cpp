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
  const std::string integralOption = "I";

  struct Promptness {
      std::string name_;
      std::string histo_name_;
  };

  std::vector<Promptness> promptnesses{
          {"prompt",    "Prompt"},
          {"nonprompt", "NonPrompt"}
  };

  TFile* fileCutVar = OpenFileWithNullptrCheck(fileNameCutVar);
  TFile* fileMC = isMc ? OpenFileWithNullptrCheck(fileNameMC) : nullptr;

  for (size_t iP = 0, nP = promptnesses.size(); iP < nP; ++iP) {
    const std::string promptness = promptnesses.at(iP).name_;

    TH1* hCutVar = GetObjectWithNullptrCheck<TH1>(fileCutVar, "hCorrYields" + promptnesses.at(iP).histo_name_);
    hCutVar->UseCurrentStyle();
    hCutVar->SetLineColor(kRed);
    hCutVar->SetMarkerColor(kRed);
    hCutVar->GetYaxis()->SetMoreLogLabels();

    TH1* hMC = isMc ? GetObjectWithNullptrCheck<TH1>(fileMC, "gen/" + promptness + "/hT") : nullptr;
    if (isMc) {
      hMC->UseCurrentStyle();
      hMC->SetLineColor(kBlue);
      hMC->SetMarkerColor(kBlue);
      hMC->GetYaxis()->SetMoreLogLabels();

      hMC = dynamic_cast<TH1D*>(hMC->Rebin(lifetimeRanges.size() - 1, hMC->GetName(),lifetimeRanges.data()));

      CustomizeHistogramsYRange({hCutVar, hMC}, true);
    } // isMc

    auto PrintYields = [&](const std::string& ccName, TH1* hCutVar, TH1* hMc, const std::string& promptness, const std::string& priBra) {
      TCanvas cc("cc", "");
      cc.SetLogy();
      cc.SetCanvasSize(1200, 800);
      if (isMc) hMc->Draw("E1");
      hCutVar->Draw("E1 same");
      const double legy1 = isMc ? 0.7 : 0.8;
      TLegend legPrompt(0.7, legy1, 0.9, 0.9);
      legPrompt.SetHeader(promptness.c_str());
      if (isMc) legPrompt.AddEntry(hMc, "MC", "PL");
      legPrompt.AddEntry(hCutVar, "Cut Var", "PL");
      legPrompt.Draw("same");
      cc.Print((ccName + ".pdf" + priBra).c_str(), "pdf");
    };
    PrintYields("compar", hCutVar, hMC, promptness, EvaluatePrintingBracket(promptnesses, iP));

    if(promptness != "prompt") continue;

    TH1* hCutVarDiff = dynamic_cast<TH1*>(hCutVar->Clone());
    TH1* hMcDiff = isMc ? dynamic_cast<TH1*>(hMC->Clone()) : nullptr;
    for(const auto& h : {hCutVarDiff, hMcDiff}) {
      if(h == nullptr) continue;
      h->Scale(1.0, "width");
      h->GetYaxis()->SetTitle("dN/dT (ps^{-1})");
    }

    PrintYields("ctfit", hCutVarDiff, hMcDiff, "", "(");

    TF1* fitCutVar = FitLifetimeHisto(hCutVarDiff, integralOption);
    TF1* fitMc = isMc ? FitLifetimeHisto(hMcDiff, integralOption) : nullptr;

    auto FitResults = [](const TF1* fitFunc, const std::string& text="") {
      const std::string lifetimeFitValue = "#tau_{#Lambda_{c}} [" + text + "] = (" +
                                           to_string_with_precision(fitFunc->GetParameter(1)*1000, 1) +
                                           " #pm " +
                                           to_string_with_precision(fitFunc->GetParError(1)*1000, 1) +
                                           ") fs";
      const std::string chi2Value = "#chi^{2} / ndf [" + text + "]= " +
                                    to_string_with_significant_figures(fitFunc->GetChisquare(), 3) +
                                    " / " +
                                    std::to_string(fitFunc->GetNDF());

      return std::make_pair(lifetimeFitValue, chi2Value);
    };

    const std::pair<std::string, std::string> fitResultsCutVar = FitResults(fitCutVar, "rec");
    const std::pair<std::string, std::string> fitResultsMC = isMc ? FitResults(fitMc, "MC") : std::pair<std::string, std::string>();
    const std::string lifetimePdg = "#tau_{#Lambda_{c}} [PDG] = (202.6 #pm 1.0) fs";
    
    const float textX1 = 0.73;
    const float textX2 = 0.86;
    const float textY2 = 0.78;
    const float textYStep = 0.06;

    TLegend leg(textX1, textY2, textX2, textY2 + 2*textYStep);
    if (isMc) leg.AddEntry(hMcDiff, "MC", "PL");
    leg.AddEntry(hCutVarDiff, "rec", "PL");

    TCanvas ccFit("ccFit", "");
    ccFit.SetCanvasSize(1200, 800);
    ccFit.SetLogy();
    if (isMc) {
      hMcDiff->Draw("E1");
      fitMc->Draw("same");
      AddOneLineText(fitResultsMC.first, {textX1, textY2 - 1*textYStep, textX2, textY2});
      AddOneLineText(fitResultsMC.second, {textX1, textY2 - 2*textYStep, textX2, textY2 - 1*textYStep});
    }
    hCutVarDiff->Draw("E1 same");
    fitCutVar->Draw("same");
    AddOneLineText(fitResultsCutVar.first, {textX1, textY2 - 3*textYStep, textX2, textY2 - 2*textYStep});
    AddOneLineText(fitResultsCutVar.second, {textX1, textY2 - 4*textYStep, textX2, textY2 - 3*textYStep});
    AddOneLineText(lifetimePdg, {textX1, textY2 - 5*textYStep, textX2, textY2 - 4*textYStep});
    leg.Draw("same");
    ccFit.Print("ctfit.pdf", "pdf");
    
    TH1* hCutVarRatio = dynamic_cast<TH1*>(hCutVarDiff->Clone());
    TH1* hMcRatio = isMc ? dynamic_cast<TH1*>(hMcDiff->Clone()) : nullptr;
    ScalePlotVertically(hCutVarRatio, hCutVarDiff, 2);
    hCutVarRatio->GetYaxis()->SetTitle("Data / Fit");
    DivideFunctionByHisto(hCutVarRatio, fitCutVar, integralOption);
    if(isMc) {
      ScalePlotVertically(hMcRatio, hMcDiff, 2);
      hMcRatio->GetYaxis()->SetTitle("Data / Fit");
      DivideFunctionByHisto(hMcRatio, fitMc, integralOption);
      CustomizeHistogramsYRange({hCutVarRatio, hMcRatio});
    }
    TCanvas ccRatio("ccRatio", "");
    ScaleCanvasVertically(&ccRatio, &ccFit, 2);
    TF1 oneline("oneline", "[0]", 0, 2);
    oneline.SetParameter(0, 1);
    oneline.SetLineColor(kBlack);
    oneline.SetLineStyle(7);
    if(isMc) hMcRatio->Draw("HIST PE");
    hCutVarRatio->Draw("HIST PE same");
    oneline.Draw("same");
    ccRatio.Print("ctfit.pdf)", "pdf");
  } // promptnesses
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./corrected_yields_qa2 fileNameCutVar (fileNameMC="")" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileNameCutVar = argv[1];
  const std::string fileNameMC = argc > 2 ? argv[2] : "";

  corrected_yields_qa2(fileNameCutVar, fileNameMC);

  return 0;
}