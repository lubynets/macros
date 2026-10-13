//
// Created by oleksii on 13.10.2026
//

#ifndef QA2_HELPERCTFIT_HPP
#define QA2_HELPERCTFIT_HPP

#include <TGraphErrors.h>
#include <TH1.h>

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

}

#endif //QA2_HELPERCTFIT_HPP
