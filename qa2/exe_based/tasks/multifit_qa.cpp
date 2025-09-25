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
#include <vector>

using namespace HelperGeneral;
using namespace HelperMath;
using namespace HelperPlot;

void MultiFitQa() {
  LoadMacro("styles/mc_qa2.style.cc");
  const std::string histoName = "hRawYieldsSignal";
  const std::string fileNameTemplate = "RawYields_Lc/RawYields_Lc";
  std::vector<double> bdtScores;
  for(int i=0; i<=2; i++) {
    bdtScores.emplace_back(0.01*i);
  }

  std::vector<int> trialNumbers{1, 2};

  const size_t nBdtScores = bdtScores.size();
  TFile* fileMarkUp = OpenFileWithNullptrCheck(std::to_string(trialNumbers.at(0)) + "/" + fileNameTemplate + ".NPgt0.01.root");
  TH1* histoMarkUp = GetObjectWithNullptrCheck<TH1>(fileMarkUp, histoName);
  const size_t nLifetimeRanges = histoMarkUp->GetNbinsX();
  fileMarkUp->Close();

  tensor<TGraphErrors*, 2> graph = make_tensor<TGraphErrors*, 2>({nLifetimeRanges, nBdtScores}, nullptr);
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
    for(size_t iScore=0; iScore<nBdtScores; ++iScore){
      const std::string fileName = std::to_string(trial) + "/" + fileNameTemplate + ".NPgt" + to_string_with_precision(bdtScores.at(iScore), 2) + ".root";
      TFile* fileIn = OpenFileWithNullptrCheck(fileName);
      TH1* histoIn = GetObjectWithNullptrCheck<TH1>(fileIn, histoName);
      for(int iBin=1; iBin<=nLifetimeRanges; ++iBin) {
        const double value = histoIn->GetBinContent(iBin);
        const double error = histoIn->GetBinError(iBin);
        auto gr = graph.at(iBin-1).at(iScore);
        gr->AddPoint(trial, value);
        gr->SetPointError(gr->GetN()-1, 0, error);
      }
      fileIn->Close();
    }
  }

  for(size_t iT=0; iT<nLifetimeRanges; ++iT) {
    for(size_t iScore=0; iScore<nBdtScores; ++iScore) {
      const std::string priBra = EvaluatePrintingBracket(nBdtScores, iScore);
      TCanvas cc("cc", "");
      cc.SetCanvasSize(1200, 800);
      graph.at(iT).at(iScore)->Draw("APE");
      cc.Print(("T_" + std::to_string(iT) + ".pdf" + priBra).c_str(), "pdf");
    }
  }
}

int main(int argc, char* argv[]) {

  MultiFitQa();
}