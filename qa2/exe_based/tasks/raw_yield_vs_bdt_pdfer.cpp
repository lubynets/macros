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
std::pair<double, double> EvaluateMeanAndStdDevOfGraph(const TGraph* graph, double from=-1e9, double to=1e9);

void raw_yield_vs_bdt_pdfer(const std::string& fileNameTemplate, const std::string& histoName) {
  LoadMacro("styles/mc_qa2.style.cc");
  //=================================================================
  const std::vector<std::string> targetSignals {/*"P" , */"NP"};
  std::vector<double> bdtScores;
  for(int i=0; i<=99; i++) {
    bdtScores.emplace_back(0.01*i);
  }
  const int moveAverageRadius = 2;
  const int moveAverageExcludeOwnPoint = true;
  const double ratioSigmaTolerance = 1.;
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
    kRatioSigma,
    nGraphs
  };

  Helper::tensor3<TGraphErrors*> graph(nGraphs);
  for(auto& gr1 : graph) {
    gr1.resize(lifeTimeRanges.size()-1);
    for(auto& gr2 : gr1) {
      gr2.resize(targetSignals.size());
    }
  }

  for(int iTargetSignal=0; iTargetSignal<targetSignals.size(); ++iTargetSignal) {
    for(int iLifeTimeRange = 0; iLifeTimeRange<lifeTimeRanges.size()-1; ++iLifeTimeRange) {
      auto InitGraph = [&](TGraphErrors*& gra, const std::string& prefixTitle, const std::string& yAxisTitle, Color_t color) {
        gra = new TGraphErrors();
        gra->SetName((prefixTitle + "_vs_" + targetSignals.at(iTargetSignal) + "_T" + std::to_string(iLifeTimeRange)).c_str());
        gra->SetMarkerColor(color);
        gra->SetLineColor(color);
        gra->SetLineWidth(1);
        gra->GetYaxis()->SetTitle(yAxisTitle.c_str());
        gra->SetTitle(("bin #" + std::to_string(iLifeTimeRange) + "#; T#in (" + to_string_with_precision(lifeTimeRanges.at(iLifeTimeRange), 2) + "#; " + to_string_with_precision(lifeTimeRanges.at(iLifeTimeRange+1), 2) + ") ps").c_str());
        gra->SetMarkerSize(0.6);
        gra->GetXaxis()->SetTitle(("bdt score " + targetSignals.at(iTargetSignal)).c_str());
        gra->GetXaxis()->SetNdivisions(315);
      };

      InitGraph(graph.at(kYield).at(iLifeTimeRange).at(iTargetSignal), "grYield", "Raw yield", kBlack);
      InitGraph(graph.at(kYieldError).at(iLifeTimeRange).at(iTargetSignal), "grYieldErr", "Raw yield error", kBlack);
      InitGraph(graph.at(kChi2).at(iLifeTimeRange).at(iTargetSignal), "grChi2", "#chi^{2}/ndf", kBlack);
      InitGraph(graph.at(kMoveAve).at(iLifeTimeRange).at(iTargetSignal), "grMoveAve", "Raw yield (move ave)", kRed);
      InitGraph(graph.at(kRatio).at(iLifeTimeRange).at(iTargetSignal), "grRatio", "Raw yield / Moving average", kBlack);
      InitGraph(graph.at(kRatioSigma).at(iLifeTimeRange).at(iTargetSignal), "grRatioSigma", "(Ratio-1)/sigma)", kBlack);
    } // lifeTimeRanges
  } // targetSignals

  for(const auto& score : bdtScores) {
    for(int iTargetSignal=0; iTargetSignal<targetSignals.size(); ++iTargetSignal) {
      TFile* fileIn = OpenFileWithNullptrCheck(fileNameTemplate + "." + targetSignals.at(iTargetSignal) + "gt" + to_string_with_precision(score, 2) + ".root");
      TH1* histoYield = GetObjectWithNullptrCheck<TH1>(fileIn, histoName);
      TH1* histoChi2 = GetObjectWithNullptrCheck<TH1>(fileIn, "hRawYieldsChiSquare");
      for(int iLifeTimeRange=0; iLifeTimeRange<lifeTimeRanges.size()-1; ++iLifeTimeRange) {
        auto gr = graph.at(kYield).at(iLifeTimeRange).at(iTargetSignal);
        gr->SetPoint(gr->GetN(), score, histoYield->GetBinContent(iLifeTimeRange + 1));
        gr->SetPointError(gr->GetN()-1, 0, histoYield->GetBinError(iLifeTimeRange + 1));

        auto gre = graph.at(kYieldError).at(iLifeTimeRange).at(iTargetSignal);
        gre->SetPoint(gre->GetN(), score, histoYield->GetBinError(iLifeTimeRange + 1));

        auto grc = graph.at(kChi2).at(iLifeTimeRange).at(iTargetSignal);
        grc->SetPoint(grc->GetN(), score, histoChi2->GetBinContent(iLifeTimeRange + 1));
      } // lifeTimeRanges
      fileIn->Close();
      ++iTargetSignal;
    } // targetSignals
  } // bdtScores

  for(int iTargetSignal=0; iTargetSignal<targetSignals.size(); ++iTargetSignal) {
    std::set<double> outlierBdts;
    for(int iLifeTimeRange=0; iLifeTimeRange<lifeTimeRanges.size()-1; ++iLifeTimeRange) {
      const std::string priBra = lifeTimeRanges.size()-1 == 1 ? "" : iLifeTimeRange == 0 ? "(" : iLifeTimeRange == lifeTimeRanges.size()-2 ? ")" : "";
      auto gry = graph.at(kYield).at(iLifeTimeRange).at(iTargetSignal);
      auto grm = graph.at(kMoveAve).at(iLifeTimeRange).at(iTargetSignal);
      auto grr = graph.at(kRatio).at(iLifeTimeRange).at(iTargetSignal);
      auto grs = graph.at(kRatioSigma).at(iLifeTimeRange).at(iTargetSignal);
      EvaluateMovingAverage(gry, grm, moveAverageRadius, moveAverageExcludeOwnPoint);
      for(int iPoint=0, nPoints = gry->GetN(); iPoint<nPoints; ++iPoint) {
        grr->SetPoint(iPoint, gry->GetPointX(iPoint), gry->GetPointY(iPoint));
      }
      DivideGraph(grr, grm);
      grr->GetYaxis()->SetRangeUser(0.9, 1.1);
      auto [meanRatio, sigmaRatio] = EvaluateMeanAndStdDevOfGraph(grr, 0.03, 0.61);
      for(int iPoint=0, nPoints = grr->GetN(); iPoint<nPoints; ++iPoint) {
        const double x = grr->GetPointX(iPoint);
        const double y = (grr->GetPointY(iPoint) - meanRatio) / sigmaRatio;
        grs->SetPoint(iPoint, x, y);
        if(std::fabs(y) > ratioSigmaTolerance) outlierBdts.insert(x);
      }
      grs->GetYaxis()->SetRangeUser(-5, 5);

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
      PrintCanvas({graph.at(kRatioSigma)}, "RatioSigma");
    } // lifeTimeRanges
    std::cout << "\nOutliers for " << targetSignals.at(iTargetSignal) << " :\n";
    for(const auto& ob : outlierBdts) {
      std::cout << ob << ", ";
    }
    std::cout << "\n\nGood bdt scores for " << targetSignals.at(iTargetSignal) << " :\n";
    for(const auto& score : bdtScores) {
      if(outlierBdts.find(score) == outlierBdts.end()) std::cout << score << ", ";
    }
    std::cout << "\n\n";
  } // targetSignals
}

std::pair<double, double> EvaluateMeanAndStdDevOfGraph(const TGraph* graph, double from, double to) {
  int nPoints{0};
  double sumValues{0.};
  double sumValues2{0.};
  for(int iPoint=0, N=graph->GetN(); iPoint<N; ++iPoint) {
    const double x = graph->GetPointX(iPoint);
    if(x<from || x>to) continue;
    const double y = graph->GetPointY(iPoint);
    sumValues += y;
    sumValues2 += y*y;
    ++nPoints;
  }
  const double mean = sumValues / nPoints;
  const double variance = sumValues2/nPoints - mean*mean;
  const double sigma = std::sqrt(variance*nPoints/(nPoints-1));

  return std::make_pair(mean, sigma);
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