//
// Created by oleksii on 13.10.2026
//

#include "HelperCtFit.hpp"
#include "HelperGeneral.hpp"

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH1.h>
#include <TLegend.h>
#include <TLine.h>

#include <array>
#include <string>
#include <utility>
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

std::pair<std::string, std::string> HelperCtFit::FitResults(const TF1* func, const std::string& text) {
  const std::string lifetimeFitValue = "#tau_{#Lambda_{c}} " + text + "= (" +
                                      HelperGeneral::to_string_with_precision(func->GetParameter(1)*1000, 1) +
                                      " #pm " +
                                      HelperGeneral::to_string_with_precision(func->GetParError(1)*1000, 1) +
                                      ") fs";
  const std::string chi2Value = "#chi^{2} / ndf = " +
                                HelperGeneral::to_string_with_significant_figures(func->GetChisquare(), 3) +
                                " / " +
                                std::to_string(func->GetNDF());

  return std::make_pair(lifetimeFitValue, chi2Value);
}

void HelperCtFit::PrintCanvasDrops(TGraphErrors* graph, int nBins) {
  const double leftEdge = -2.;
  const double rightEdge = graph->GetPointX(graph->GetN()-1) + 2;

  TCanvas ccDrops("ccDrops", "");
  ccDrops.SetCanvasSize(1200, 800);
  graph->GetXaxis()->SetLimits(leftEdge, rightEdge);
  graph->Draw("AP");
  TGraphErrors* grMean = GetSubGraph(graph, SubGraphType::MeanPoint, kRed, nBins);
  TGraphErrors* grOneDrop = GetSubGraph(graph, SubGraphType::OneDrop, kBlue, nBins);
  TGraphErrors* grTwoFit = GetSubGraph(graph, SubGraphType::TwoFit, kGreen+2, nBins);
  grMean->Draw("PE same");
  grOneDrop->Draw("PE same");
  grTwoFit->Draw("PE same");
  TLine lineMean(leftEdge, graph->GetPointY(0), rightEdge, graph->GetPointY(0));
  TLine lineUp(leftEdge, graph->GetPointY(0)+graph->GetErrorY(0), rightEdge, graph->GetPointY(0)+graph->GetErrorY(0));
  TLine lineDown(leftEdge, graph->GetPointY(0)-graph->GetErrorY(0), rightEdge, graph->GetPointY(0)-graph->GetErrorY(0));
  lineUp.SetLineStyle(7);
  lineDown.SetLineStyle(7);
  for(const auto& line : {&lineMean, &lineUp, &lineDown}) {
    line->SetLineWidth(2);
    line->SetLineColor(kRed);
    line->Draw("same");
  }
  TLegend leg(0.75, 0.70, 0.95, 0.90);
  leg.AddEntry(grMean, "Fit all points", "PE");
  leg.AddEntry(grOneDrop, "Fit w/o one point", "PE");
  leg.AddEntry(grTwoFit, "Fit two points", "PE");
  leg.Draw("same");

  ccDrops.Print("dropSummary.pdf", "pdf");
}

std::array<TGraphErrors*, static_cast<int>(HelperCtFit::GraphQaType::NGraphQas)> HelperCtFit::CreateQaGraphs() {
  std::array<TGraphErrors*, static_cast<int>(GraphQaType::NGraphQas)> graphs;
  for(auto& gr : graphs) {
    gr = new TGraphErrors();
    gr->GetXaxis()->SetTitle("drop set");
  }
  graphs.at(static_cast<int>(GraphQaType::Tau))->SetName("grTau");
  graphs.at(static_cast<int>(GraphQaType::Integral))->SetName("grA");
  graphs.at(static_cast<int>(GraphQaType::Chi2))->SetName("grChi2");
  graphs.at(static_cast<int>(GraphQaType::Eta))->SetName("grEta");
  graphs.at(static_cast<int>(GraphQaType::Tau))->GetYaxis()->SetTitle("#tau (fs)");
  graphs.at(static_cast<int>(GraphQaType::Integral))->GetYaxis()->SetTitle("A");
  graphs.at(static_cast<int>(GraphQaType::Chi2))->GetYaxis()->SetTitle("#chi^{2} / ndf");
  graphs.at(static_cast<int>(GraphQaType::Eta))->GetYaxis()->SetTitle("#eta, %");

  return graphs;
}

void HelperCtFit::FillQaGraphs(std::array<TGraphErrors*, static_cast<int>(GraphQaType::NGraphQas)>& graphs, const TF1* func, int dropSet, double referenceSigmaTau) {
  const auto binsToDrop = EvalBinsToDrop(dropSet);

  graphs.at(static_cast<int>(GraphQaType::Integral))->AddPoint(dropSet, func->GetParameter(0)*1000);
  graphs.at(static_cast<int>(GraphQaType::Integral))->SetPointError(graphs.at(static_cast<int>(GraphQaType::Integral))->GetN()-1, 0., func->GetParError(0)*1000);
  graphs.at(static_cast<int>(GraphQaType::Tau))->AddPoint(dropSet, func->GetParameter(1)*1000);
  graphs.at(static_cast<int>(GraphQaType::Tau))->SetPointError(graphs.at(static_cast<int>(GraphQaType::Tau))->GetN()-1, 0., func->GetParError(1)*1000);

  const auto chi2 = func->GetChisquare();
  const auto ndf = func->GetNDF();
  if(ndf != 0) graphs.at(static_cast<int>(GraphQaType::Chi2))->AddPoint(dropSet, chi2 / ndf);

  if (dropSet == 0) referenceSigmaTau = func->GetParError(1);
  if(binsToDrop.size() < 2) {
    const double currentSigmaTau = func->GetParError(1);
    const double eta = (currentSigmaTau / referenceSigmaTau - 1.) * 100.;
    graphs.at(static_cast<int>(GraphQaType::Eta))->AddPoint(dropSet, eta);
  }
}
