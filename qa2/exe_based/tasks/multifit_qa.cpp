//
// Created by oleksii on 25.09.2025.
//
#include "HelperGeneral.hpp"
#include "HelperMath.hpp"
#include "HelperPlot.hpp"

#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>

#include <iostream>
#include <numeric>
#include <vector>

using namespace HelperGeneral;
using namespace HelperMath;
using namespace HelperPlot;

void MultiFitQa() {
  LoadMacro("styles/mc_qa2.style.cc");
  const bool isVerbose{false};
  const std::string histoName = "hRawYieldsSignal";
  const std::string fileNameTemplate = "RawYields_Lc/RawYields_Lc";
  const int nTrials = 100;
  std::vector<double> bdtScores;
  for(int i=0; i<=99; i++) {
    bdtScores.emplace_back(0.01*i);
  }

  std::vector<int> trialNumbers(nTrials);
  std::iota(trialNumbers.begin(), trialNumbers.end(), 1);

  const size_t nBdtScores = bdtScores.size();
  TFile* fileMarkUp = OpenFileWithNullptrCheck(std::to_string(trialNumbers.at(0)) + "/" + fileNameTemplate + ".NPgt0.01.root");
  TH1* histoMarkUp = GetObjectWithNullptrCheck<TH1>(fileMarkUp, histoName);
  const size_t nLifetimeRanges = histoMarkUp->GetNbinsX();

  tensor<TGraphErrors*, 2> graph = make_tensor<TGraphErrors*, 2>({nLifetimeRanges, nBdtScores}, nullptr);
  tensor<std::map<double, size_t>, 2> values = make_tensor<std::map<double, size_t>, 2>({nLifetimeRanges, nBdtScores}, {});
  tensor<std::map<double, size_t>, 2> errors = make_tensor<std::map<double, size_t>, 2>({nLifetimeRanges, nBdtScores}, {});
  for(size_t iT=0; iT<nLifetimeRanges; ++iT) {
    for(size_t iScore=0; iScore<nBdtScores; ++iScore) {
      auto& gr = graph.at(iT).at(iScore);
      gr = new TGraphErrors();
      gr->SetName(("gr_T" + std::to_string(iT) + "_NPgt" + to_string_with_precision(bdtScores.at(iScore), 2)).c_str());
      gr->SetTitle(("gr_T" + std::to_string(iT) + "_NPgt" + to_string_with_precision(bdtScores.at(iScore), 2)).c_str());
      gr->GetXaxis()->SetTitle("Trial #");
      gr->GetYaxis()->SetTitle(histoName.c_str());
    }
  }

  for(const auto& trial : trialNumbers) {
    for(size_t iScore=0; iScore<nBdtScores; ++iScore) {
      const std::string fileName = std::to_string(trial) + "/" + fileNameTemplate + ".NPgt" + to_string_with_precision(bdtScores.at(iScore), 2) + ".root";
      if(isVerbose) std::cout << "Opening " << fileName;
      TFile* fileIn = TFile::Open(fileName.c_str(), "read");
      if(isVerbose) std::cout << ", opened successfully\n";
      if(fileIn == nullptr) {
        std::cout << fileName << " is missing\n";
        continue;
      }
      if(isVerbose) std::cout << "Reading histogram";
      TH1* histoIn = GetObjectWithNullptrCheck<TH1>(fileIn, histoName);
      if(isVerbose) std::cout << ", read successfully\n";
      if(isVerbose) std::cout << "iBin = ";
      for(int iBin=1; iBin<=nLifetimeRanges; ++iBin) {
        if(isVerbose) std::cout << iBin << " ";
        const double value = histoIn->GetBinContent(iBin);
        const double error = histoIn->GetBinError(iBin);
        auto gr = graph.at(iBin-1).at(iScore);
        gr->AddPoint(trial, value);
        gr->SetPointError(gr->GetN()-1, 0, error);
        values.at(iBin-1).at(iScore).insert({value, trial});
        errors.at(iBin-1).at(iScore).insert({error, trial});
      } // nLifetimeRanges
      if(isVerbose) std::cout << "\nColsing " << fileName;
      fileIn->Close();
      if(isVerbose) std::cout << ", closed successfully\n\n";
    } // nBdtScores
  } // trialNumbers

  auto FindMapMedian = [](const std::map<double, size_t>& map) {
    const int mapSize = map.size();
    auto it = map.begin();
    std::advance(it, mapSize/2);
    return it->second;
  };

  auto FindGraphsPointByX = [](const TGraph* gr, double x) {
    const int nPoints = gr->GetN();
    for(int iPoint=0; iPoint<nPoints; ++ iPoint) {
      const double grX = gr->GetPointX(iPoint);
      if(std::fabs(grX - x) < 1e-4) return iPoint;
    }
    return -1;
  };

  for(size_t iScore=0; iScore<nBdtScores; ++iScore) {
    const std::string priBra = EvaluatePrintingBracket(nBdtScores, iScore);
    TH1* histoSmooth = dynamic_cast<TH1*>(histoMarkUp->Clone());
    histoSmooth->Reset();
    for(size_t iT=0; iT<nLifetimeRanges; ++iT) {
      TCanvas cc("cc", "");
      cc.SetCanvasSize(1200, 800);
      const auto& gr = graph.at(iT).at(iScore);
      gr->Draw("APE");
      const int medianValueTrial = FindMapMedian(values.at(iT).at(iScore));
      const int medianValuePoint = FindGraphsPointByX(gr, medianValueTrial);
      const double medianValue = gr->GetPointY(medianValuePoint);
      const int medianErrorTrial = FindMapMedian(errors.at(iT).at(iScore));
      const int medianErrorPoint = FindGraphsPointByX(gr, medianErrorTrial);
      const double medianError = gr->GetErrorY(medianErrorPoint);
      histoSmooth->SetBinContent(iT+1, medianValue);
      histoSmooth->SetBinError(iT+1, medianError);
      TF1* lineValue = HorizontalLine4Graph(medianValue, gr);
      TF1* lineErrorUp = HorizontalLine4Graph(medianValue+medianError, gr);
      TF1* lineErrorDown = HorizontalLine4Graph(medianValue-medianError, gr);
      lineErrorUp->SetLineStyle(7);
      lineErrorDown->SetLineStyle(7);
      for(const auto& line : {lineValue, lineErrorUp, lineErrorDown}) {
        line->SetLineColor(kRed);
        line->SetLineWidth(2);
        line->Draw("same");
      }
      cc.Print(("T_" + std::to_string(iT) + ".pdf" + priBra).c_str(), "pdf");
    } // nLifetimeRanges
    TFile* fileSmooth = TFile::Open(("smooth/" + fileNameTemplate + ".NPgt" + to_string_with_precision(bdtScores.at(iScore), 2) + ".root").c_str(), "recreate");
    histoSmooth->Write();
    fileSmooth->Close();
  } // nBdtScores
  fileMarkUp->Close();
}

int main(int argc, char* argv[]) {

  MultiFitQa();
}