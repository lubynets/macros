//
// Created by oleksii on 22.06.25.
//
#include "Helper.hpp"

#include <TFile.h>
#include <TGraphErrors.h>

#include <iostream>
#include <vector>

using namespace Helper;

void raw_yield_vs_bdt_pdfer(const std::string& fileNameTemplate, const std::string& histoName) {
  LoadMacro("styles/mc_qa2.style.cc");
  //=================================================================
  const std::vector<std::string> targetSignals {"P" , "NP"};
  std::vector<float> bdtScores;
  for(int i=0; i<=99; i++) {
    bdtScores.emplace_back(0.01*i);
  }
  //=================================================================

  TFile* fileMarkup = OpenFileWithNullptrCheck(fileNameTemplate + ".Pgt" + to_string_with_precision(bdtScores.at(0), 2) + ".root");
  TH1* histoMarkup = GetObjectWithNullptrCheck<TH1>(fileMarkup, histoName);
  std::vector<float> lifeTimeRanges;
  for(int iBin=1; iBin<=histoMarkup->GetNbinsX()+1; ++iBin) {
    lifeTimeRanges.emplace_back(histoMarkup->GetBinLowEdge(iBin));
  }
  fileMarkup->Close();

  std::vector<std::vector<TGraphErrors*>> grYield(lifeTimeRanges.size()-1);
  for(auto& gr : grYield) {
    gr.resize(targetSignals.size());
  }
  std::vector<std::vector<TGraph*>> grYieldErr(lifeTimeRanges.size()-1);
  for(auto& gr : grYieldErr) {
    gr.resize(targetSignals.size());
  }

  int iTargetSignal{0};
  for(const auto& tarSig : targetSignals) {
    for(int iLifeTimeRange = 0; iLifeTimeRange<lifeTimeRanges.size()-1; ++iLifeTimeRange) {
      grYield.at(iLifeTimeRange).at(iTargetSignal) = new TGraphErrors(bdtScores.size());
      auto gr = grYield.at(iLifeTimeRange).at(iTargetSignal);
      gr->SetName(("grYield_vs_" + tarSig + "_T" + std::to_string(iLifeTimeRange)).c_str());
      gr->SetTitle(("bin #" + std::to_string(iLifeTimeRange) + "#; T#in (" + to_string_with_precision(lifeTimeRanges.at(iLifeTimeRange), 2) + "#; " + to_string_with_precision(lifeTimeRanges.at(iLifeTimeRange+1), 2) + ") ps").c_str());
      gr->SetMarkerColor(kBlack);
      gr->SetLineColor(kBlack);
      gr->SetMarkerSize(0.6);
      gr->SetLineWidth(1);
      gr->GetXaxis()->SetTitle(("bdt score " + tarSig).c_str());
      gr->GetYaxis()->SetTitle("Raw yield");
      gr->GetXaxis()->SetNdivisions(315);

      grYieldErr.at(iLifeTimeRange).at(iTargetSignal) = new TGraph(bdtScores.size());
      auto gre = grYieldErr.at(iLifeTimeRange).at(iTargetSignal);
      gre->SetName(("grYieldErr_vs_" + tarSig + "_T" + std::to_string(iLifeTimeRange)).c_str());
      gre->SetTitle(("bin #" + std::to_string(iLifeTimeRange) + "#; T#in (" + to_string_with_precision(lifeTimeRanges.at(iLifeTimeRange), 2) + "#; " + to_string_with_precision(lifeTimeRanges.at(iLifeTimeRange+1), 2) + ") ps").c_str());
      gre->SetMarkerColor(kBlack);
      gre->SetMarkerSize(0.6);
      gre->GetXaxis()->SetTitle(("bdt score " + tarSig).c_str());
      gre->GetYaxis()->SetTitle("Raw yield error");
      gre->GetXaxis()->SetNdivisions(315);
    } // lifeTimeRanges
    ++iTargetSignal;
  } // targetSignals

  for(const auto& score : bdtScores) {
    iTargetSignal = 0;
    for(const auto& tarSig: targetSignals) {
      TFile* fileIn = OpenFileWithNullptrCheck(fileNameTemplate + "." + tarSig + "gt" + to_string_with_precision(score, 2) + ".root");
      TH1* histoIn = GetObjectWithNullptrCheck<TH1>(fileIn, histoName);
      for(int iLifeTimeRange=0; iLifeTimeRange<lifeTimeRanges.size()-1; ++iLifeTimeRange) {
        auto gr = grYield.at(iLifeTimeRange).at(iTargetSignal);
        gr->SetPoint(gr->GetN(), score, histoIn->GetBinContent(iLifeTimeRange+1));
        gr->SetPointError(gr->GetN()-1, 0, histoIn->GetBinError(iLifeTimeRange+1));
        auto gre = grYieldErr.at(iLifeTimeRange).at(iTargetSignal);
        gre->SetPoint(gr->GetN(), score, histoIn->GetBinError(iLifeTimeRange+1));
      }
      fileIn->Close();
      ++iTargetSignal;
    } // targetSignals
  } // bdtScores

  for(int iLifeTimeRange=0; iLifeTimeRange<lifeTimeRanges.size()-1; ++iLifeTimeRange) {
    const std::string priBra = lifeTimeRanges.size()-1 == 1 ? "" : iLifeTimeRange == 0 ? "(" : iLifeTimeRange == lifeTimeRanges.size()-2 ? ")" : "";
    for(int iTargetSignal=0; iTargetSignal<targetSignals.size(); ++iTargetSignal) {
      TCanvas ccYield("ccYield", "");
      ccYield.SetGridx();
      ccYield.SetCanvasSize(1200, 800);
      grYield.at(iLifeTimeRange).at(iTargetSignal)->Draw("APE");
      ccYield.Print(("grRawYield_vs_" + targetSignals.at(iTargetSignal) + "." + histoName + ".pdf" + priBra).c_str());

      TCanvas ccYieldErr("ccYieldErr", "");
      ccYieldErr.SetGridx();
      ccYieldErr.SetCanvasSize(1200, 800);
      grYieldErr.at(iLifeTimeRange).at(iTargetSignal)->Draw("AP");
      ccYieldErr.Print(("grErr_vs_" + targetSignals.at(iTargetSignal) + "." + histoName + ".pdf" + priBra).c_str());
    } // targetSignals
  } // lifeTimeRanges
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./raw_yield_vs_bdt_pdfer fileNameTemplate (histoName=hRawYields)" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileNameTemplate = argv[1];
  const std::string histoName = argc > 2 ? argv[2] : "hRawYields";

  raw_yield_vs_bdt_pdfer(fileNameTemplate, histoName);

  return 0;
}