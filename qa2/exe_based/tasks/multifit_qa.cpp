//
// Created by oleksii on 25.09.2025.
//
#include "HelperGeneral.hpp"
#include "HelperMath.hpp"
#include "HelperPlot.hpp"

#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TStyle.h>

#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>
#include <unordered_set>
#include <utility>

using namespace HelperGeneral;
using namespace HelperMath;
using namespace HelperPlot;

constexpr bool IsVerbose{false};

void MultiFitQa(const std::string& strategy) {
  if(strategy != "median" && strategy != "chi2" && strategy != "medianSmart") throw std::runtime_error("MultiFitQa(): strategy must be 'median', 'chi2' or 'medianSmart'");

  LoadMacro("styles/mc_qa2.style.cc");
  gStyle->SetMarkerSize(1);
  const std::string fileNameTemplate = "RawYields_Lc/RawYields_Lc";
  const int nTrials = 200;
  std::vector<double> bdtScores;
  for(int i=1; i<=99; i++) {
    bdtScores.emplace_back(0.01*i);
  }
//   bdtScores={0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90};

  std::vector<std::string> variables {
   "hRawYieldsSignal",
   "hRawYieldsSignalCounted",
   "hRawYieldsSigma",
   "hRawYieldsMean",
//    "hRawYieldsSignificance",
//    "hRawYieldsSgnOverBkg",
//    "hRawYieldsBkg",
//    "hRawYieldsChiSquareBkg",
//    "hRawYieldsDscbAlphaL",
//    "hRawYieldsDscbAlphaR",
//    "hRawYieldsDscbNL",
//    "hRawYieldsDscbNR"
   "hRawYieldsChiSquareTotal"
  };

  std::vector<int> trialNumbers(nTrials);
  std::iota(trialNumbers.begin(), trialNumbers.end(), 1);

  const size_t nBdtScores = bdtScores.size();

  TFile* fileMarkUp{nullptr};
  for(const auto& trial : trialNumbers) {
    fileMarkUp = TFile::Open(("trials/" + std::to_string(trial) + "/" + fileNameTemplate + ".NPgt" + to_string_with_precision(bdtScores.at(0), 2) + ".root").c_str());
    if(fileMarkUp != nullptr) break;
  }
  if(fileMarkUp == nullptr) throw std::runtime_error("fileMarkUp == nullptr");

  TH1* histoMarkUp{nullptr};
  int iExistingVar{0};
  size_t nLifetimeRanges{0};
  for(int iVar=0, nVars=variables.size(); iVar<nVars; ++iVar) {
    histoMarkUp = fileMarkUp->Get<TH1>(variables.at(iExistingVar).c_str());
    if(histoMarkUp == nullptr) {
      variables.erase(variables.begin() + iExistingVar);
      continue;
    }
    if(nLifetimeRanges == 0) nLifetimeRanges = histoMarkUp->GetNbinsX();
    ++iExistingVar;
  }

  const size_t nVars = variables.size();

  const auto itChi2 = std::find(variables.begin(), variables.end(), "hRawYieldsChiSquareTotal");
  const auto iChi2 = std::distance(variables.begin(), itChi2);
  if((strategy == "chi2" || strategy == "medianSmart") && iChi2 == nVars) throw std::runtime_error("MultiFitQa(): strategy is '" + strategy + "', but hRawYieldsChiSquareTotal is missing");

  tensor<TGraphErrors*, 3> graph = make_tensor<TGraphErrors*, 3>({nVars, nLifetimeRanges, nBdtScores}, nullptr);
  tensor<std::map<double, size_t>, 3> values = make_tensor<std::map<double, size_t>, 3>({nVars, nLifetimeRanges, nBdtScores}, {});
  tensor<std::map<double, size_t>, 3> errors = make_tensor<std::map<double, size_t>, 3>({nVars, nLifetimeRanges, nBdtScores}, {});
  tensor<TGraphErrors*, 3> graphVsChi2 = make_tensor<TGraphErrors*, 3>({nVars, nLifetimeRanges, nBdtScores}, nullptr);
  tensor<TGraphErrors*, 3> graphErr = make_tensor<TGraphErrors*, 3>({nVars, nLifetimeRanges, nBdtScores}, nullptr);
  tensor<TGraphErrors*, 3> graphErrVsChi2 = make_tensor<TGraphErrors*, 3>({nVars, nLifetimeRanges, nBdtScores}, nullptr);
  for(size_t iVar=0; iVar<nVars; ++iVar) {
    for (size_t iT = 0; iT < nLifetimeRanges; ++iT) {
      for (size_t iScore = 0; iScore < nBdtScores; ++iScore) {
        auto& gr = graph.at(iVar).at(iT).at(iScore);
        gr = new TGraphErrors();
        const std::string grName = "T" + std::to_string(iT) + "_NPgt" + to_string_with_precision(bdtScores.at(iScore), 2);
        gr->SetName((variables.at(iVar) + "_" + grName).c_str());
        gr->SetTitle(grName.c_str());
        gr->GetXaxis()->SetTitle("Trial #");
        gr->GetYaxis()->SetTitle(variables.at(iVar).c_str());

        auto& grVsChi2 = graphVsChi2.at(iVar).at(iT).at(iScore);
        grVsChi2 = new TGraphErrors();
        grVsChi2->SetName((variables.at(iVar) + "_" + grName).c_str());
        grVsChi2->SetTitle(grName.c_str());
        grVsChi2->GetXaxis()->SetTitle("#chi^{2}/ndf");
        grVsChi2->GetYaxis()->SetTitle(variables.at(iVar).c_str());

        auto& grErr = graphErr.at(iVar).at(iT).at(iScore);
        grErr = new TGraphErrors();
        const std::string grErrName = "T" + std::to_string(iT) + "_NPgt" + to_string_with_precision(bdtScores.at(iScore), 2);
        grErr->SetName((variables.at(iVar) + "_" + grName).c_str());
        grErr->SetTitle(grName.c_str());
        grErr->GetXaxis()->SetTitle("Trial #");
        grErr->GetYaxis()->SetTitle((variables.at(iVar) + " error").c_str());

        auto& grErrVsChi2 = graphErrVsChi2.at(iVar).at(iT).at(iScore);
        grErrVsChi2 = new TGraphErrors();
        grErrVsChi2->SetName((variables.at(iVar) + "_" + grName).c_str());
        grErrVsChi2->SetTitle(grName.c_str());
        grErrVsChi2->GetXaxis()->SetTitle("#chi^{2}/ndf");
        grErrVsChi2->GetYaxis()->SetTitle((variables.at(iVar) + " error").c_str());
      } // nBdtScores
    } // nLifetimeRanges
  } // nVars

  for(const auto& trial : trialNumbers) {
    for(size_t iScore=0; iScore<nBdtScores; ++iScore) {
      const std::string fileName = "trials/" + std::to_string(trial) + "/" + fileNameTemplate + ".NPgt" + to_string_with_precision(bdtScores.at(iScore), 2) + ".root";
      if(IsVerbose) std::cout << "Opening " << fileName;
      TFile* fileIn = TFile::Open(fileName.c_str(), "read");
      if(IsVerbose) std::cout << ", opened successfully\n";
      if(fileIn == nullptr) {
        if(IsVerbose) std::cout << ", " << fileName << " is missing\n";
        continue;
      }
      for(size_t iVar=0; iVar<nVars; ++iVar) {
        if(IsVerbose) std::cout << "Reading histogram" << variables.at(iVar);
        TH1* histoIn = fileIn->Get<TH1>(variables.at(iVar).c_str());
        TH1* histoChi2 = fileIn->Get<TH1>("hRawYieldsChiSquareTotal");
        if(histoIn == nullptr || histoChi2 == nullptr) {
          if(IsVerbose) std::cout << ", " << variables.at(iVar) << " or hRawYieldsChiSquareTotal is missing\n";
          continue;
        }
        if(IsVerbose) std::cout << ", read successfully\t";
        if(IsVerbose) std::cout << "iBin = ";
        for (int iBin = 1; iBin <= static_cast<int>(nLifetimeRanges); ++iBin) {
          if(IsVerbose) std::cout << iBin << " ";
          const double value = histoIn->GetBinContent(iBin);
          const double error = histoIn->GetBinError(iBin);
          auto gr = graph.at(iVar).at(iBin - 1).at(iScore);
          gr->AddPoint(trial, value);
          gr->SetPointError(gr->GetN() - 1, 0, error);
          graphErr.at(iVar).at(iBin - 1).at(iScore)->AddPoint(trial, error);
          values.at(iVar).at(iBin - 1).at(iScore).insert({value, trial});
          errors.at(iVar).at(iBin - 1).at(iScore).insert({error, trial});

          const double chi2 = histoChi2->GetBinContent(iBin);
          auto grVsChi2 = graphVsChi2.at(iVar).at(iBin - 1).at(iScore);
          grVsChi2->AddPoint(chi2, value);
          grVsChi2->SetPointError(grVsChi2->GetN() - 1, 0, error);
          graphErrVsChi2.at(iVar).at(iBin - 1).at(iScore)->AddPoint(chi2, error);
        } // nLifetimeRanges
        if(IsVerbose) std::cout << "\n";
      } // nVars
      if(IsVerbose) std::cout << "\nClosing " << fileName;
      fileIn->Close();
      if(IsVerbose) std::cout << ", closed successfully\n\n";
    } // nBdtScores
  } // trialNumbers

  auto FindMapMedian = [](const std::map<double, size_t>& map) {
    const size_t mapSize = map.size();
    auto it = map.begin();
    std::advance(it, mapSize/2);
    return it->second;
  };

  auto FindBadTrials = [](const std::map<double, size_t>& map) {
    const size_t mapSize = map.size();
    auto it = map.begin();
    std::advance(it, mapSize/5);
    const double criticalChi2 = it->first;

    std::vector<int> badTrials{};
    for (; it != map.end(); ++it) {
      badTrials.push_back(it->second);
    }

    return std::pair<std::vector<int>, double>{badTrials, criticalChi2};
  };

  auto RemoveBadTrials = [](std::map<double, size_t>& map, const std::vector<int> badTrials) {
    std::unordered_set<int> badTrialsSet(badTrials.begin(), badTrials.end());

    for(auto it = map.begin(); it != map.end(); ) {
      if(badTrialsSet.count(it->second)) it = map.erase(it);
      else ++it;
    }
  };

  auto FindGraphsPointByX = [](const TGraph* gr, double x) {
    const int nPoints = gr->GetN();
    for(int iPoint=0; iPoint<nPoints; ++ iPoint) {
      const double grX = gr->GetPointX(iPoint);
      if(std::fabs(grX - x) < 1e-4) return iPoint;
    }
    return -1;
  };

  for (size_t iScore = 0; iScore < nBdtScores; ++iScore) {
    const std::string priBra = EvaluatePrintingBracket(nBdtScores, iScore);
    const std::string fileSmoothName = "smooth/" + fileNameTemplate + ".NPgt" + to_string_with_precision(bdtScores.at(iScore), 2) + ".root";
    auto getDirectory = [](const std::string &s) -> std::string {
      const auto pos = s.rfind('/');  // find last slash
      if (pos == std::string::npos) {
        return "";  // no slash found -> return empty string
      }
      return s.substr(0, pos);
    };
    const std::string dirSmoothPath = getDirectory(fileSmoothName);
    if(!dirSmoothPath.empty()) MkDirBash(dirSmoothPath);
    TFile* fileSmooth = TFile::Open(("smooth/" + fileNameTemplate + ".NPgt" + to_string_with_precision(bdtScores.at(iScore), 2) + ".root").c_str(), "recreate");
    for(size_t iVar=0; iVar<nVars; ++iVar) {
      TH1* histoSmooth = dynamic_cast<TH1*>(GetObjectWithNullptrCheck<TH1>(fileMarkUp, variables.at(iVar))->Clone());
      histoSmooth->Reset();
      for (size_t iT = 0; iT < nLifetimeRanges; ++iT) {
        TCanvas cc("cc", "", 1200, 800);
        TCanvas ccErr("ccErr", "", 1200, 800);
        cc.SetGrid();
        ccErr.SetGrid();
        TCanvas ccVsChi2("ccVsChi2", "", 1200, 800);
        TCanvas ccErrVsChi2("ccErrVsChi2", "", 1200, 800);
        const auto& gr = graph.at(iVar).at(iT).at(iScore);
        const auto& grVsChi2 = graphVsChi2.at(iVar).at(iT).at(iScore);
        cc.cd();
        gr->Draw("APE");
        ccVsChi2.cd();
        grVsChi2->Draw("APE");
        ccErr.cd();
        graphErr.at(iVar).at(iT).at(iScore)->Draw("APE");
        ccErrVsChi2.cd();
        graphErrVsChi2.at(iVar).at(iT).at(iScore)->Draw("APE");

        const auto [badTrials, criticalChi2] = strategy == "medianSmart" ? FindBadTrials(values.at(iChi2).at(iT).at(iScore)) : std::pair<std::vector<int>, double>{std::vector<int>{}, UndefValueDouble};
        if(iVar != iChi2) {
          RemoveBadTrials(values.at(iVar).at(iT).at(iScore), badTrials);
          RemoveBadTrials(errors.at(iVar).at(iT).at(iScore), badTrials);
        }

        const size_t bestValueTrial = strategy == "median" || strategy == "medianSmart" ? FindMapMedian(values.at(iVar).at(iT).at(iScore)) : values.at(iChi2).at(iT).at(iScore).begin()->second;
        const int bestValuePoint = FindGraphsPointByX(gr, static_cast<double>(bestValueTrial));
        const double bestValue = gr->GetPointY(bestValuePoint);
        const size_t bestErrorTrial = strategy == "median" || strategy == "medianSmart" ? FindMapMedian(errors.at(iVar).at(iT).at(iScore)) : values.at(iChi2).at(iT).at(iScore).begin()->second;
        const int bestErrorPoint = FindGraphsPointByX(gr, static_cast<double>(bestErrorTrial));
        const double bestError = gr->GetErrorY(bestErrorPoint);
        histoSmooth->SetBinContent(iT + 1, bestValue);
        histoSmooth->SetBinError(iT + 1, bestError);
        gr->SetTitle((static_cast<std::string>(gr->GetTitle()) + " (" + std::to_string(bestValueTrial) + ", " + std::to_string(bestErrorTrial) + ")").c_str());
        grVsChi2->SetTitle((static_cast<std::string>(grVsChi2->GetTitle()) + " (" + std::to_string(bestValueTrial) + ", " + std::to_string(bestErrorTrial) + ")" + ", chi2 crit. = " + to_string_with_precision(criticalChi2, 1)).c_str());
        graphErr.at(iVar).at(iT).at(iScore)->SetTitle((static_cast<std::string>(graphErr.at(iVar).at(iT).at(iScore)->GetTitle()) + " (" + std::to_string(bestValueTrial) + ", " + std::to_string(bestErrorTrial) + ")").c_str());
        graphErrVsChi2.at(iVar).at(iT).at(iScore)->SetTitle((static_cast<std::string>(graphErrVsChi2.at(iVar).at(iT).at(iScore)->GetTitle()) + " (" + std::to_string(bestValueTrial) + ", " + std::to_string(bestErrorTrial) + ")" + ", chi2 crit. = " + to_string_with_precision(criticalChi2, 1)).c_str());
        TF1* lineValue = HorizontalLine4Graph(bestValue, gr);
        TF1* lineErrorUp = HorizontalLine4Graph(bestValue + bestError, gr);
        TF1* lineErrorDown = HorizontalLine4Graph(bestValue - bestError, gr);
        TF1* lineValueChi2 = HorizontalLine4Graph(bestValue, grVsChi2);
        TF1* lineErrorUpChi2 = HorizontalLine4Graph(bestValue + bestError, grVsChi2);
        TF1* lineErrorDownChi2 = HorizontalLine4Graph(bestValue - bestError, grVsChi2);
        TF1* lineError = HorizontalLine4Graph(bestError, gr);
        TF1* lineErrorChi2 = HorizontalLine4Graph(bestError, grVsChi2);
        for (const auto& line: {lineErrorUp, lineErrorDown, lineErrorUpChi2, lineErrorDownChi2}) {
          line->SetLineStyle(7);
          line->SetLineStyle(7);
        }
        int iLine{0};
        for (const auto& line: {lineValue, lineErrorUp, lineErrorDown, lineValueChi2, lineErrorUpChi2, lineErrorDownChi2, lineError, lineErrorChi2}) {
          line->SetLineColor(kRed);
          line->SetLineWidth(2);
          if(iLine<3) cc.cd();
          if(iLine >=3 && iLine < 6) ccVsChi2.cd();
          if(iLine == 6) ccErr.cd();
          if(iLine == 7) ccErrVsChi2.cd();
          line->Draw("same");
          ++iLine;
        }
        MkDirBash("hTrials/" + variables.at(iVar));
        MkDirBash("hTrials/" + variables.at(iVar) + "VsChi2");
        MkDirBash("hTrials/err." + variables.at(iVar));
        MkDirBash("hTrials/err." + variables.at(iVar) + "VsChi2");
        cc.Print(("hTrials/" + variables.at(iVar) + "/" + variables.at(iVar) + "_T_" + std::to_string(iT+1) + ".pdf" + priBra).c_str(), "pdf");
        ccVsChi2.Print(("hTrials/" + variables.at(iVar) + "VsChi2" + "/" + variables.at(iVar) + "VsChi2" + "_T_" + std::to_string(iT+1) + ".pdf" + priBra).c_str(), "pdf");
        ccErr.Print(("hTrials/err." + variables.at(iVar) + "/" + variables.at(iVar) + "_T_" + std::to_string(iT+1) + ".pdf" + priBra).c_str(), "pdf");
        ccErrVsChi2.Print(("hTrials/err." + variables.at(iVar) + "VsChi2" + "/" + variables.at(iVar) + "VsChi2" + "_T_" + std::to_string(iT+1) + ".pdf" + priBra).c_str(), "pdf");
      } // nLifetimeRanges
      fileSmooth->cd();
      histoSmooth->Write();
    } // nVars
     fileSmooth->Close();
  } // nBdtScores
  fileMarkUp->Close();
}

int main(int argc, char* argv[]) {

  const std::string strategy = argc > 1 ? argv[1] : "median";

  MultiFitQa(strategy);
}
