//
// Created by oleksii on 13.10.2026
//

#include "HelperCtFit.hpp"

#include <TGraphErrors.h>
#include <TH1.h>

#include <vector>

void HelperCtFit::ExcludeBin(TH1* h, int binNumber) {
  h->SetBinContent(binNumber, 0.);
  h->SetBinError(binNumber, 0.);
}

std::vector<int> HelperCtFit::EvalBinsToDrop(int dropSet) {
  std::vector<int> result{};
  int iBit{0};
  while(dropSet >> iBit > 0) {
    if((dropSet >> iBit) & 1) result.push_back(iBit+1);
    ++iBit;
  }

  return result;
}

int HelperCtFit::EvalDropSet(const std::vector<int>& bins) {
  int result = 0;

  for(const auto& bin : bins) {
    result |= (1 << (bin - 1));
  }

  return result;
}

TGraphErrors* HelperCtFit::GetSubGraph(TGraphErrors* grIn, SubGraphType subGraphType, Color_t color, int nFittedPoints) {
  TGraphErrors* grOut = new TGraphErrors();
  grOut->SetName(grIn->GetName());
  for(int iPoint=0, nPoints=grIn->GetN(); iPoint<nPoints; ++iPoint) {
    const int dropSet = static_cast<int>(std::round(grIn->GetPointX(iPoint)));
    if(subGraphType == SubGraphType::MeanPoint && dropSet !=0) continue;
    if(subGraphType == SubGraphType::OneDrop && !((dropSet > 0) && ((dropSet & (dropSet - 1)) == 0))) continue; // dropSet is a power of 2, i.e. 1, 2, 4, 8, ...
    if(subGraphType == SubGraphType::TwoFit && __builtin_popcount(dropSet) != nFittedPoints-2) continue;
    grOut->AddPoint(grIn->GetPointX(iPoint), grIn->GetPointY(iPoint));
    grOut->SetPointError(grOut->GetN()-1, 0., grIn->GetErrorY(iPoint));
  }
  grOut->SetMarkerColor(color);
  grOut->SetLineColor(color);

  return grOut;
}
