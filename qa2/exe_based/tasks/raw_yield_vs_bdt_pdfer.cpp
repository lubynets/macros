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
  const std::vector<std::string> targetSignals {/*"P" , */"NP"};
  std::vector<float> bdtScores;
  for(int i=0; i<=99; i++) {
    bdtScores.emplace_back(0.01*i);
  }
  //=================================================================

  TFile* fileMarkup = OpenFileWithNullptrCheck(fileNameTemplate + "." + targetSignals.at(0) + "gt" + to_string_with_precision(bdtScores.at(0), 2) + ".root");
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
  std::vector<std::vector<TGraphErrors*>> grYieldErr(lifeTimeRanges.size()-1);
  for(auto& gr : grYieldErr) {
    gr.resize(targetSignals.size());
  }
  std::vector<std::vector<TGraphErrors*>> grChi2(lifeTimeRanges.size()-1);
  for(auto& gr : grChi2) {
    gr.resize(targetSignals.size());
  }

  int iTargetSignal{0};
  for(const auto& tarSig : targetSignals) {
    for(int iLifeTimeRange = 0; iLifeTimeRange<lifeTimeRanges.size()-1; ++iLifeTimeRange) {
      grYield.at(iLifeTimeRange).at(iTargetSignal) = new TGraphErrors(bdtScores.size());
      auto gr = grYield.at(iLifeTimeRange).at(iTargetSignal);
      gr->SetName(("grYield_vs_" + tarSig + "_T" + std::to_string(iLifeTimeRange)).c_str());
      gr->SetMarkerColor(kBlack);
      gr->SetLineColor(kBlack);
      gr->SetLineWidth(1);
      gr->GetYaxis()->SetTitle("Raw yield");

      grYieldErr.at(iLifeTimeRange).at(iTargetSignal) = new TGraphErrors(bdtScores.size());
      auto gre = grYieldErr.at(iLifeTimeRange).at(iTargetSignal);
      gre->SetName(("grYieldErr_vs_" + tarSig + "_T" + std::to_string(iLifeTimeRange)).c_str());
      gre->SetMarkerColor(kBlack);
      gre->GetYaxis()->SetTitle("Raw yield error");

      grChi2.at(iLifeTimeRange).at(iTargetSignal) = new TGraphErrors(bdtScores.size());
      auto grc = grChi2.at(iLifeTimeRange).at(iTargetSignal);
      grc->SetName(("grChi2_vs_" + tarSig + "_T" + std::to_string(iLifeTimeRange)).c_str());
      grc->SetMarkerColor(kBlack);
      grc->GetYaxis()->SetTitle("#chi^{2}/ndf");

      for(const auto& g : std::vector<TGraphErrors*>{gr, gre, grc}) {
        g->SetTitle(("bin #" + std::to_string(iLifeTimeRange) + "#; T#in (" + to_string_with_precision(lifeTimeRanges.at(iLifeTimeRange), 2) + "#; " + to_string_with_precision(lifeTimeRanges.at(iLifeTimeRange+1), 2) + ") ps").c_str());
        g->SetMarkerSize(0.6);
        g->GetXaxis()->SetTitle(("bdt score " + tarSig).c_str());
        g->GetXaxis()->SetNdivisions(315);
      }
    } // lifeTimeRanges
    ++iTargetSignal;
  } // targetSignals

  for(const auto& score : bdtScores) {
    iTargetSignal = 0;
    for(const auto& tarSig: targetSignals) {
      TFile* fileIn = OpenFileWithNullptrCheck(fileNameTemplate + "." + tarSig + "gt" + to_string_with_precision(score, 2) + ".root");
      TH1* histoYield = GetObjectWithNullptrCheck<TH1>(fileIn, histoName);
      TH1* histoChi2 = GetObjectWithNullptrCheck<TH1>(fileIn, "hRawYieldsChiSquare");
      for(int iLifeTimeRange=0; iLifeTimeRange<lifeTimeRanges.size()-1; ++iLifeTimeRange) {
        auto gr = grYield.at(iLifeTimeRange).at(iTargetSignal);
        gr->SetPoint(gr->GetN(), score, histoYield->GetBinContent(iLifeTimeRange + 1));
        gr->SetPointError(gr->GetN()-1, 0, histoYield->GetBinError(iLifeTimeRange + 1));

        auto gre = grYieldErr.at(iLifeTimeRange).at(iTargetSignal);
        gre->SetPoint(gre->GetN(), score, histoYield->GetBinError(iLifeTimeRange + 1));

        auto grc = grChi2.at(iLifeTimeRange).at(iTargetSignal);
        grc->SetPoint(grc->GetN(), score, histoChi2->GetBinContent(iLifeTimeRange + 1));
      }
      fileIn->Close();
      ++iTargetSignal;
    } // targetSignals
  } // bdtScores

  for(int iLifeTimeRange=0; iLifeTimeRange<lifeTimeRanges.size()-1; ++iLifeTimeRange) {
    const std::string priBra = lifeTimeRanges.size()-1 == 1 ? "" : iLifeTimeRange == 0 ? "(" : iLifeTimeRange == lifeTimeRanges.size()-2 ? ")" : "";
    for(int iTargetSignal=0; iTargetSignal<targetSignals.size(); ++iTargetSignal) {
      auto PrintCanvas = [&](const std::vector<std::vector<TGraphErrors*>>& gr, const std::string& name) {
        TCanvas cc(("cc" + name).c_str(), "");
        cc.SetGridx();
        cc.SetCanvasSize(1200, 800);
        gr.at(iLifeTimeRange).at(iTargetSignal)->Draw("APE");
        cc.Print(("gr" + name + "_vs_" + targetSignals.at(iTargetSignal) + "." + histoName + ".pdf" + priBra).c_str());
      };

      PrintCanvas(grYield, "RawYield");
      PrintCanvas(grYieldErr, "Err");
      PrintCanvas(grChi2, "Chi2");
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