//
// Created by oleksii on 22.06.25.
//
#include "Helper.hpp"

#include <TFile.h>
#include <TGraphErrors.h>
#include <TStyle.h>

#include <iostream>
#include <set>
#include <vector>

using namespace Helper;
std::pair<double, double> EvaluateMeanAndStdDevOfGraph(const TGraph* graph, double from=-1e9, double to=1e9);

void raw_yield_vs_bdt_pdfer(const std::string& fileNameTemplate, const std::string& histoName) {
  LoadMacro("styles/mc_qa2.style.cc");
  gStyle->SetMarkerSize(0.6);
  gStyle->SetNdivisions(315, "X");
  gStyle->SetLineWidth(1);
  //=================================================================
  const std::vector<std::string> targetSignals {/*"P" , */"NP"};
  std::vector<double> bdtScores;
  for(int i=0; i<=99; i++) {
    bdtScores.emplace_back(0.01*i);
  }
  const int moveAverageRadius = 3;
  const int moveAverageExcludeOwnPoint = true;
  const double ratioSigmaTolerance = 2.;
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
    nGraphs
  };

  enum AdvancedGraph : short {
    kMoveAve = 0,
    kRatio,
    kRatioSigma,
    nAdvancedGraphs
  };

  Helper::tensor3<TGraphErrors*> graph = Helper::CreateTensor3<TGraphErrors*>(nGraphs, lifeTimeRanges.size()-1, targetSignals.size());

  for(int iTargetSignal=0; iTargetSignal<targetSignals.size(); ++iTargetSignal) {
    for(int iLifeTimeRange = 0; iLifeTimeRange<lifeTimeRanges.size()-1; ++iLifeTimeRange) {
      auto InitGraph = [&](TGraphErrors*& gra, const std::string& prefixTitle, const std::string& yAxisTitle) {
        gra = new TGraphErrors();
        gra->SetName((prefixTitle + "_vs_" + targetSignals.at(iTargetSignal) + "_T" + std::to_string(iLifeTimeRange)).c_str());
        gra->GetYaxis()->SetTitle(yAxisTitle.c_str());
        gra->SetTitle(("bin #" + std::to_string(iLifeTimeRange) + "#; T#in (" + to_string_with_precision(lifeTimeRanges.at(iLifeTimeRange), 2) + "#; " + to_string_with_precision(lifeTimeRanges.at(iLifeTimeRange+1), 2) + ") ps").c_str());
        gra->GetXaxis()->SetTitle(("bdt score " + targetSignals.at(iTargetSignal)).c_str());
      };

      InitGraph(graph.at(kYield).at(iLifeTimeRange).at(iTargetSignal), "grYield", "Raw yield");
      InitGraph(graph.at(kYieldError).at(iLifeTimeRange).at(iTargetSignal), "grYieldErr", "Raw yield error");
      InitGraph(graph.at(kChi2).at(iLifeTimeRange).at(iTargetSignal), "grChi2", "#chi^{2}/ndf");
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
    auto ProcessMoveAverageRatioSigma = [&](TGraphErrors*& grIn, std::array<TGraphErrors*, 3>& graph) {
      EvaluateMovingAverage(grIn, graph.at(kMoveAve), moveAverageRadius, moveAverageExcludeOwnPoint);
      for(int iPoint=0, nPoints = grIn->GetN(); iPoint<nPoints; ++iPoint) {
        graph.at(kRatio)->SetPoint(iPoint, grIn->GetPointX(iPoint), grIn->GetPointY(iPoint));
      }
      DivideGraph(graph.at(kRatio), graph.at(kMoveAve));
      auto [meanRatio, sigmaRatio] = EvaluateMeanAndStdDevOfGraph(graph.at(kRatio), 0.03, 0.61);
      for(int iPoint=0, nPoints = graph.at(kRatio)->GetN(); iPoint<nPoints; ++iPoint) {
        const double x = graph.at(kRatio)->GetPointX(iPoint);
        const double y = (graph.at(kRatio)->GetPointY(iPoint) - meanRatio) / sigmaRatio;
        graph.at(kRatioSigma)->SetPoint(iPoint, x, y);
        if(std::fabs(y) > ratioSigmaTolerance) outlierBdts.insert(x);
      }
      graph.at(kRatio)->GetYaxis()->SetRangeUser(0.9, 1.1);
      graph.at(kRatioSigma)->GetYaxis()->SetRangeUser(-5, 5);
    };
    auto InitMoveAverageRatioSigma = [&](const std::string& baseName) {
      TGraphErrors* grm = new TGraphErrors();
      TGraphErrors* grr = new TGraphErrors();
      TGraphErrors* grs = new TGraphErrors();
      grm->SetMarkerColor(kRed);
      for(const auto& gr : {grm, grr, grs}) {
        gr->GetXaxis()->SetTitle(("bdt score " + targetSignals.at(iTargetSignal)).c_str());
      }
      grm->GetYaxis()->SetTitle(("Moving average (" + baseName + ")").c_str());
      grr->GetYaxis()->SetTitle((baseName + " / Moving average").c_str());
      grs->GetYaxis()->SetTitle(("(Ratio-1)/sigma (" + baseName + ")").c_str());

      return std::array<TGraphErrors*, 3>{grm, grr, grs};
    };
    auto PrintCanvas = [&](const std::vector<TGraphErrors*>& gr, const std::string& name, const std::string& priBra) {
      TCanvas cc(("cc" + name).c_str(), "");
      cc.SetGridx();
      cc.SetCanvasSize(1200, 800);
      int iGr{0};
      for(const auto& grr : gr) {
        const std::string option = iGr == 0 ? "APE" : "PE same";
        grr->Draw(option.c_str());
        ++iGr;
      }
      cc.Print(("gr" + name + "_vs_" + targetSignals.at(iTargetSignal) + "." + histoName + ".pdf" + priBra).c_str());
    };

    for(int iLifeTimeRange=0; iLifeTimeRange<lifeTimeRanges.size()-1; ++iLifeTimeRange) {
      const std::string priBra = lifeTimeRanges.size()-1 == 1 ? "" : iLifeTimeRange == 0 ? "(" : iLifeTimeRange == lifeTimeRanges.size()-2 ? ")" : "";

      auto gry = graph.at(kYield).at(iLifeTimeRange).at(iTargetSignal);
      auto gre = graph.at(kYieldError).at(iLifeTimeRange).at(iTargetSignal);

      auto grYieldAdv = InitMoveAverageRatioSigma("Raw yield");
      ProcessMoveAverageRatioSigma(gry, grYieldAdv);

      auto grErrorAdv = InitMoveAverageRatioSigma("Raw yield error");
      ProcessMoveAverageRatioSigma(gre, grErrorAdv);
      
      PrintCanvas({gry, grYieldAdv.at(kMoveAve)}, "RawYield", priBra);
      PrintCanvas({grYieldAdv.at(kRatio)}, "RawYieldRatio", priBra);
      PrintCanvas({grYieldAdv.at(kRatioSigma)}, "RawYieldRatioSigma", priBra);

      PrintCanvas({gre, grErrorAdv.at(kMoveAve)}, "Err", priBra);
      PrintCanvas({grErrorAdv.at(kRatio)}, "ErrRatio", priBra);
      PrintCanvas({grErrorAdv.at(kRatioSigma)}, "ErrRatioSigma", priBra);

      PrintCanvas({graph.at(kChi2).at(iLifeTimeRange).at(iTargetSignal)}, "Chi2", priBra);

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
    std::cout << " ./raw_yield_vs_bdt_pdfer fileNameTemplate (histoName=hRawYieldsSignal)" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileNameTemplate = argv[1];
  const std::string histoName = argc > 2 ? argv[2] : "hRawYieldsSignal";

  raw_yield_vs_bdt_pdfer(fileNameTemplate, histoName);

  return 0;
}