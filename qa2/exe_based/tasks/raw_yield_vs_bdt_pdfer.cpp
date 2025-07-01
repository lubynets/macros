//
// Created by oleksii on 22.06.25.
//
#include "Helper.hpp"

#include <TFile.h>
#include <TGraphErrors.h>

#include <iostream>
#include <set>
#include <vector>

using namespace Helper;

void raw_yield_vs_bdt_pdfer(const std::string& fileNameTemplate, const std::string& histoName) {
  LoadMacro("styles/mc_qa2.style.cc");
  //=================================================================
  const std::vector<std::string> targetSignals {/*"P" , */"NP"};
  std::vector<double> bdtScores;
  for(int i=0; i<=99; i++) {
    bdtScores.emplace_back(0.01*i);
  }
  const int moveAverageLength = 5;
  const int moveAverageExcludeOwnPoint = true;
  const double ratioTolerance = 0.03;
  //=================================================================

  TFile* fileMarkup = OpenFileWithNullptrCheck(fileNameTemplate + "." + targetSignals.at(0) + "gt" + to_string_with_precision(bdtScores.at(0), 2) + ".root");
  TH1* histoMarkup = GetObjectWithNullptrCheck<TH1>(fileMarkup, histoName);
  std::vector<float> lifeTimeRanges;
  for(int iBin=1; iBin<=histoMarkup->GetNbinsX()+1; ++iBin) {
    lifeTimeRanges.emplace_back(histoMarkup->GetBinLowEdge(iBin));
  }
  fileMarkup->Close();

  enum Graph : short {
    kYield = 0,
    kYieldError,
    kChi2,
    kMoveAve,
    kRatio,
    nGraphs
  };

  Helper::tensor3<TGraphErrors*> graph(nGraphs);
  for(auto& gr1 : graph) {
    gr1.resize(lifeTimeRanges.size()-1);
    for(auto& gr2 : gr1) {
      gr2.resize(targetSignals.size());
    }
  }

  int iTargetSignal{0};
  for(const auto& tarSig : targetSignals) {
    for(int iLifeTimeRange = 0; iLifeTimeRange<lifeTimeRanges.size()-1; ++iLifeTimeRange) {
      auto InitGraph = [&](TGraphErrors*& gra, const std::string& prefixTitle, const std::string& yAxisTitle, Color_t color) {
        gra = new TGraphErrors();
        gra->SetName((prefixTitle + "_vs_" + tarSig + "_T" + std::to_string(iLifeTimeRange)).c_str());
        gra->SetMarkerColor(color);
        gra->SetLineColor(color);
        gra->SetLineWidth(1);
        gra->GetYaxis()->SetTitle(yAxisTitle.c_str());
        gra->SetTitle(("bin #" + std::to_string(iLifeTimeRange) + "#; T#in (" + to_string_with_precision(lifeTimeRanges.at(iLifeTimeRange), 2) + "#; " + to_string_with_precision(lifeTimeRanges.at(iLifeTimeRange+1), 2) + ") ps").c_str());
        gra->SetMarkerSize(0.6);
        gra->GetXaxis()->SetTitle(("bdt score " + tarSig).c_str());
        gra->GetXaxis()->SetNdivisions(315);
      };

      InitGraph(graph.at(kYield).at(iLifeTimeRange).at(iTargetSignal), "grYield", "Raw yield", kBlack);
      InitGraph(graph.at(kYieldError).at(iLifeTimeRange).at(iTargetSignal), "grYieldErr", "Raw yield error", kBlack);
      InitGraph(graph.at(kChi2).at(iLifeTimeRange).at(iTargetSignal), "grChi2", "#chi^{2}/ndf", kBlack);
      InitGraph(graph.at(kMoveAve).at(iLifeTimeRange).at(iTargetSignal), "grMoveAve", "Raw yield (move ave)", kRed);
      InitGraph(graph.at(kRatio).at(iLifeTimeRange).at(iTargetSignal), "grRatio", "Raw yield / Moving average", kBlack);
    } // lifeTimeRanges
    ++iTargetSignal;
  } // targetSignals

  for(const auto& score : bdtScores) {
    iTargetSignal = 0;
    for(const auto& tarSig: targetSignals) {
      TFile* fileIn = OpenFileWithNullptrCheck(fileNameTemplate + "." + tarSig + "gt" + to_string_with_precision(score, 2) + ".root");
      TH1* histoYield = GetObjectWithNullptrCheck<TH1>(fileIn, histoName);
      TH1* histoChi2 = GetObjectWithNullptrCheck<TH1>(fileIn, "hRawYieldsChiSquare");
      std::set<double> outlierBdts; // TODO fill
      for(int iLifeTimeRange=0; iLifeTimeRange<lifeTimeRanges.size()-1; ++iLifeTimeRange) {
        auto gr = graph.at(kYield).at(iLifeTimeRange).at(iTargetSignal);
        gr->SetPoint(gr->GetN(), score, histoYield->GetBinContent(iLifeTimeRange + 1));
        gr->SetPointError(gr->GetN()-1, 0, histoYield->GetBinError(iLifeTimeRange + 1));

        auto gre = graph.at(kYieldError).at(iLifeTimeRange).at(iTargetSignal);
        gre->SetPoint(gre->GetN(), score, histoYield->GetBinError(iLifeTimeRange + 1));

        auto grc = graph.at(kChi2).at(iLifeTimeRange).at(iTargetSignal);
        grc->SetPoint(grc->GetN(), score, histoChi2->GetBinContent(iLifeTimeRange + 1));
      } // lifeTimeRanges
      for(int iLifeTimeRange=0; iLifeTimeRange<lifeTimeRanges.size()-1; ++iLifeTimeRange) {
        EvaluateMovingAverage(graph.at(kYield).at(iLifeTimeRange).at(iTargetSignal), graph.at(kMoveAve).at(iLifeTimeRange).at(iTargetSignal), moveAverageLength, moveAverageExcludeOwnPoint);
        for(int iPoint=0, nPoints = graph.at(kYield).at(iLifeTimeRange).at(iTargetSignal)->GetN(); iPoint<nPoints; ++iPoint) {
          graph.at(kRatio).at(iLifeTimeRange).at(iTargetSignal)->SetPoint(iPoint, graph.at(kYield).at(iLifeTimeRange).at(iTargetSignal)->GetPointX(iPoint), graph.at(kYield).at(iLifeTimeRange).at(iTargetSignal)->GetPointY(iPoint));
        }
        DivideGraph(graph.at(kRatio).at(iLifeTimeRange).at(iTargetSignal), graph.at(kMoveAve).at(iLifeTimeRange).at(iTargetSignal));
        graph.at(kRatio).at(iLifeTimeRange).at(iTargetSignal)->GetYaxis()->SetRangeUser(0.9, 1.1);
      }
      fileIn->Close();
      ++iTargetSignal;
    } // targetSignals
  } // bdtScores

  for(int iLifeTimeRange=0; iLifeTimeRange<lifeTimeRanges.size()-1; ++iLifeTimeRange) {
    const std::string priBra = lifeTimeRanges.size()-1 == 1 ? "" : iLifeTimeRange == 0 ? "(" : iLifeTimeRange == lifeTimeRanges.size()-2 ? ")" : "";
    for(int iTargetSignal=0; iTargetSignal<targetSignals.size(); ++iTargetSignal) {
      auto PrintCanvas = [&](const Helper::tensor3<TGraphErrors*>& gr, const std::string& name) {
        TCanvas cc(("cc" + name).c_str(), "");
        cc.SetGridx();
        cc.SetCanvasSize(1200, 800);
        int iGr{0};
        for(const auto& grr : gr) {
          const std::string option = iGr == 0 ? "APE" : "PE same";
          grr.at(iLifeTimeRange).at(iTargetSignal)->Draw(option.c_str());
          ++iGr;
        }
        cc.Print(("gr" + name + "_vs_" + targetSignals.at(iTargetSignal) + "." + histoName + ".pdf" + priBra).c_str());
      };

      PrintCanvas({graph.at(kYield), graph.at(kMoveAve)}, "RawYield");
      PrintCanvas({graph.at(kYieldError)}, "Err");
      PrintCanvas({graph.at(kChi2)}, "Chi2");
      PrintCanvas({graph.at(kRatio)}, "Ratio");
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