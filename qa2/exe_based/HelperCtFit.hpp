//
// Created by oleksii on 13.10.2026
//

#ifndef QA2_HELPERCTFIT_HPP
#define QA2_HELPERCTFIT_HPP

#include <TGraphErrors.h>
#include <TF1.h>
#include <TH1.h>

#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace HelperCtFit
{
enum class RunModes {
  MeanFitOnly = 0,
  AllWoOne,
  AllPossible
};

enum class SubGraphType {
  MeanPoint = 0,
  OneDrop,
  TwoFit
};

void ExcludeBin(TH1* h, int binNumber);

std::vector<int> EvalBinsToDrop(int dropSet);

int EvalDropSet(const std::vector<int>& bins);

TGraphErrors* GetSubGraph(TGraphErrors* grIn, SubGraphType subGraphType, Color_t color, int nFittedPoints);

std::pair<std::string, std::string> FitResults(const TF1* fitFunc, const std::string& text="");

void PrintCanvasDrops(TGraphErrors* graph, int nBins, const std::string& priBra);

constexpr std::string_view LifetimePdgSv = "#tau_{#Lambda_{c}} [PDG] = (202.6 #pm 1.0) fs";
}

#endif //QA2_HELPERCTFIT_HPP
