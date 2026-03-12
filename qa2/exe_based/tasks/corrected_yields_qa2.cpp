//
// Created by oleksii on 26.06.25.
//
#include "HelperGeneral.hpp"
#include "HelperMath.hpp"
#include "HelperPlot.hpp"

#include <TGraphErrors.h>
#include <TFile.h>
#include <TLegend.h>
#include <TLine.h>

#include <cmath>
#include <iostream>
#include <string>

using namespace HelperGeneral;
using namespace HelperMath;
using namespace HelperPlot;

constexpr bool IsSaveCanvasAsRoot{true};
constexpr bool IsDoBinsDropping{true};

void ExcludeBin(TH1* h, int binNumber);
std::vector<int> EvalBinsToDrop(int dropSet);
TGraphErrors* GetSubGraph(TGraphErrors* grIn, int subGraphType, Color_t color, int nFittedPoints);

enum SubGraphType {
  MeanPoint = 0,
  OneDrop,
  TwoFit
};

void corrected_yields_qa2(const std::string& fileNameCutVar, const std::string& fileNameMC) {
  LoadMacro("styles/mc_qa2.style.cc");

  const bool isMc = !fileNameMC.empty();

  const std::vector<double> lifetimeRanges = {0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8};
  const std::string integralOption = "I";

  struct Promptness {
    std::string name_;
    std::string histo_name_;
  };

  const std::vector<Promptness> promptnesses {
    {"prompt",    "Prompt"},
    {"nonprompt", "NonPrompt"}
  };

  TFile* fileCutVar = OpenFileWithNullptrCheck(fileNameCutVar);
  TFile* fileMC = isMc ? OpenFileWithNullptrCheck(fileNameMC) : nullptr;

  for (size_t iP = 0, nP = promptnesses.size(); iP < nP; ++iP) {
    const std::string promptness = promptnesses.at(iP).name_;

    TH1* hCutVar = GetObjectWithNullptrCheck<TH1>(fileCutVar, "hCorrYields" + promptnesses.at(iP).histo_name_);
    CheckTAxisForRanges(*hCutVar->GetXaxis(), lifetimeRanges);
    hCutVar = dynamic_cast<TH1D*>(hCutVar->Rebin(lifetimeRanges.size() - 1, hCutVar->GetName(),lifetimeRanges.data()));
    hCutVar->UseCurrentStyle();
    hCutVar->SetLineColor(kRed);
    hCutVar->SetMarkerColor(kRed);

    TH1* hMC = isMc ? GetObjectWithNullptrCheck<TH1>(fileMC, "gen/" + promptness + "/hT") : nullptr;
    if (isMc) {
      CheckTAxisForRanges(*hMC->GetXaxis(), lifetimeRanges);
      hMC->UseCurrentStyle();
      hMC->SetLineColor(kBlue);
      hMC->SetMarkerColor(kBlue);

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
    if(isMc) PrintYields("compar", hCutVar, hMC, promptness, EvaluatePrintingBracket(promptnesses, iP));

    if(promptness != "prompt") continue;

    TH1* hCutVarDiff = dynamic_cast<TH1*>(hCutVar->Clone());
    TH1* hMcDiff = isMc ? dynamic_cast<TH1*>(hMC->Clone()) : nullptr;
    for(const auto& h : {hCutVarDiff, hMcDiff}) {
      if(h == nullptr) continue;
      h->Scale(1.0, "width");
      h->GetXaxis()->SetTitle("#it{t} (ps)");
      h->GetYaxis()->SetTitle("d#it{N}/d#it{t} (ps^{-1})");
    }

    TGraphErrors grTau;
    TGraphErrors grA;

    grTau.SetName("grTau");
    grTau.GetXaxis()->SetTitle("drop set");
    grTau.GetYaxis()->SetTitle("#tau (fs)");
    grA.SetName("grA");
    grA.GetXaxis()->SetTitle("drop set");
    grA.GetYaxis()->SetTitle("A");

    const int nBins = hCutVarDiff->GetNbinsX();
    const int nDropSets = IsDoBinsDropping ? 1 << nBins : 1;

    TCanvas emptycanvas("emptycanvas", "", 1200, 800);
    emptycanvas.Print("ctfit.pdf[", "pdf");

    auto FitLifetimeHistos = [&](TH1* histoRec, TH1* histoMc, int dropSet) {
      TF1* fitCutVar = FitLifetimeHisto(histoRec, integralOption);
      TF1* fitMc = isMc ? FitLifetimeHisto(histoMc, integralOption) : nullptr;

      grA.AddPoint(dropSet, fitCutVar->GetParameter(0)*1000);
      grA.SetPointError(grA.GetN()-1, 0., fitCutVar->GetParError(0)*1000);
      grTau.AddPoint(dropSet, fitCutVar->GetParameter(1)*1000);
      grTau.SetPointError(grTau.GetN()-1, 0., fitCutVar->GetParError(1)*1000);

      auto FitResults = [](const TF1* fitFunc, const std::string& text="") {
        const std::string lifetimeFitValue = "#tau_{#Lambda_{c}} [" + text + "] = (" +
                                            to_string_with_precision(fitFunc->GetParameter(1)*1000, 1) +
                                            " #pm " +
                                            to_string_with_precision(fitFunc->GetParError(1)*1000, 1) +
                                            ") fs";
        const std::string chi2Value = "#chi^{2} / ndf = " +
                                      to_string_with_significant_figures(fitFunc->GetChisquare(), 3) +
                                      " / " +
                                      std::to_string(fitFunc->GetNDF());

        return std::make_pair(lifetimeFitValue, chi2Value);
      };

      const std::pair<std::string, std::string> fitResultsCutVar = FitResults(fitCutVar, "rec");
      const std::pair<std::string, std::string> fitResultsMC = isMc ? FitResults(fitMc, "MC") : std::pair<std::string, std::string>();
      const std::string lifetimePdg = "#tau_{#Lambda_{c}} [PDG] = (202.6 #pm 1.0) fs";

      const float textX1 = 0.70;
      const float textX2 = 0.86;
      const float textY2 = 0.86;
      const float textYStep = 0.06;

      TCanvas ccFit("ccFit", "");
      ccFit.SetCanvasSize(1200, 800);
      ccFit.SetLogy();
      if (isMc) {
        histoMc->Draw("E1");
        fitMc->Draw("same");
        AddOneLineText(fitResultsMC.first, {textX1, textY2 - 1*textYStep, textX2, textY2}, "brNDC", 0.04);
        AddOneLineText(fitResultsMC.second, {textX1, textY2 - 2*textYStep, textX2, textY2 - 1*textYStep}, "brNDC", 0.04);
      }
      histoRec->Draw("E1 same");
      const auto binsToDrop = EvalBinsToDrop(dropSet);
      std::string droppedBins{};
      for(const auto& binToDrop : binsToDrop) {
        droppedBins.append(std::to_string(binToDrop));
        droppedBins.append(", ");
      }
      if(dropSet != 0) {
        droppedBins.pop_back();
        droppedBins.pop_back();
        histoRec->SetTitle(("drop set = " + std::to_string(dropSet) + ", dropped bins: {" + droppedBins + "}").c_str());
      }
      fitCutVar->Draw("same");
      AddOneLineText(fitResultsCutVar.first, {textX1, textY2 - 3*textYStep, textX2, textY2 - 2*textYStep}, "brNDC", 0.04);
      AddOneLineText(fitResultsCutVar.second, {textX1, textY2 - 4*textYStep, textX2, textY2 - 3*textYStep}, "brNDC", 0.04);
      AddOneLineText(lifetimePdg, {textX1, textY2 - 6*textYStep, textX2, textY2 - 5*textYStep}, "brNDC", 0.04);
      const std::string& priBraFit = dropSet == 0 ? "(" : "";
      std::cout << "priBraFit = " << priBraFit << "\n";
      ccFit.Print("ctfit.pdf", "pdf");

      TFile* fileOut{nullptr};
      if(IsSaveCanvasAsRoot && dropSet == 0) {
        fileOut = TFile::Open("ctfit.root", "recreate");
        fileOut->cd();
        ccFit.Write();
      }

      TH1* hCutVarRatio = dynamic_cast<TH1*>(histoRec->Clone());
      TH1* hMcRatio = isMc ? dynamic_cast<TH1*>(histoMc->Clone()) : nullptr;
      ScalePlotVertically(hCutVarRatio, histoRec, 2);
      hCutVarRatio->GetYaxis()->SetTitle("Data / Fit");
      DivideHistoByFunction(hCutVarRatio, fitCutVar, integralOption);
      if(isMc) {
        ScalePlotVertically(hMcRatio, histoMc, 2);
        hMcRatio->GetYaxis()->SetTitle("Data / Fit");
        DivideHistoByFunction(hMcRatio, fitMc, integralOption);
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
      const std::string& priBraRatio = dropSet == nDropSets-1 ? ")" : "";
      std::cout << "priBraRatio = " << priBraRatio << "\n";
      ccRatio.Print("ctfit.pdf", "pdf");
      if(IsSaveCanvasAsRoot && dropSet == 0) {
        fileOut->cd();
        ccRatio.Write();
        fileOut->Close();
      }
    };

    for(int iDropSet=0; iDropSet<nDropSets; ++iDropSet) {
      const auto binsToDrop = EvalBinsToDrop(iDropSet);
      if(static_cast<int>(binsToDrop.size()) >= nBins-1) continue;
      std::cout << iDropSet << "\t{";
      for(const auto& v : binsToDrop) {
        std::cout << v << " ";
      }
      std::cout << "}\n";

      TH1* hReco = dynamic_cast<TH1*>(hCutVarDiff->Clone());
      TH1* hMc = isMc ? dynamic_cast<TH1*>(hMcDiff->Clone()) : nullptr;

      for(const auto& binToDrop : binsToDrop) {
        ExcludeBin(hReco, binToDrop);
        if(isMc) ExcludeBin(hMc, binToDrop);
      }

      FitLifetimeHistos(hReco, hMc, iDropSet);
    } // nDropSets
    emptycanvas.Print("ctfit.pdf]", "pdf");

    const double leftEdge = -2.;
    const double rightEdge = grTau.GetPointX(grTau.GetN()-1) + 2;

    auto PrintCanvasDrops = [&](TGraphErrors* graph, const std::string& priBra) {
      TCanvas ccDrops("ccDrops", "");
      ccDrops.SetCanvasSize(1200, 800);
      graph->GetXaxis()->SetLimits(leftEdge, rightEdge);
      graph->Draw("AP");
      TGraphErrors* grMean = GetSubGraph(graph, MeanPoint, kRed, nBins);
      TGraphErrors* grOneDrop = GetSubGraph(graph, OneDrop, kBlue, nBins);
      TGraphErrors* grTwoFit = GetSubGraph(graph, TwoFit, kGreen+2, nBins);
      grMean->Draw("PE same");
      grOneDrop->Draw("PE same");
      grTwoFit->Draw("PE same");
      TLine lineMean(leftEdge, graph->GetPointY(0), rightEdge, graph->GetPointY(0));
      TLine lineUp(leftEdge, graph->GetPointY(0)+graph->GetErrorY(0), rightEdge, graph->GetPointY(0)+graph->GetErrorY(0));
      TLine lineDown(leftEdge, graph->GetPointY(0)-graph->GetErrorY(0), rightEdge, graph->GetPointY(0)-graph->GetErrorY(0));
      lineUp.SetLineStyle(7);
      lineDown.SetLineStyle(7);
      for(const auto& line : {&lineMean, &lineUp, &lineDown}) {
        line->SetLineWidth(2);
        line->SetLineColor(kRed);
        line->Draw("same");
      }
      TLegend leg(0.75, 0.70, 0.95, 0.90);
      leg.AddEntry(grMean, "Fit all points", "PE");
      leg.AddEntry(grOneDrop, "Fit w/o one point", "PE");
      leg.AddEntry(grTwoFit, "Fit two points", "PE");
      leg.Draw("same");

      ccDrops.Print(("dropSummary.pdf" + priBra).c_str(), "pdf");
    };

    PrintCanvasDrops(&grTau, "(");
    PrintCanvasDrops(&grA, ")");

  } // promptnesses
}

void ExcludeBin(TH1* h, int binNumber) {
  h->SetBinContent(binNumber, 0.);
  h->SetBinError(binNumber, 0.);
}

std::vector<int> EvalBinsToDrop(int dropSet) {
  std::vector<int> result{};
  int iBit{0};
  while(dropSet >> iBit > 0) {
    if((dropSet >> iBit) & 1) result.push_back(iBit+1);
    ++iBit;
  }

  return result;
}

TGraphErrors* GetSubGraph(TGraphErrors* grIn, int subGraphType, Color_t color, int nFittedPoints) {
  TGraphErrors* grOut = new TGraphErrors();
  grOut->SetName(grIn->GetName());
  for(int iPoint=0, nPoints=grIn->GetN(); iPoint<nPoints; ++iPoint) {
    const int dropSet = static_cast<int>(std::round(grIn->GetPointX(iPoint)));
    if(subGraphType == MeanPoint && dropSet !=0) continue;
    if(subGraphType == OneDrop && !((dropSet > 0) && ((dropSet & (dropSet - 1)) == 0))) continue; // dropSet is a power of 2, i.e. 1, 2, 4, 8, ...
    if(subGraphType == TwoFit && __builtin_popcount(dropSet) != nFittedPoints-2) continue;
    grOut->AddPoint(grIn->GetPointX(iPoint), grIn->GetPointY(iPoint));
    grOut->SetPointError(grOut->GetN()-1, 0., grIn->GetErrorY(iPoint));
  }
  grOut->SetMarkerColor(color);
  grOut->SetLineColor(color);

  return grOut;
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
