//
// Created by oleksii on 10.02.25.
//
#include "HelperGeneral.hpp"
#include "HelperPlot.hpp"

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TROOT.h>

#include <iostream>

using namespace HelperGeneral;
using namespace HelperPlot;

void graph_qa3(const std::string& fileName1, const std::string& fileName2, const std::string& promptness, const std::string& constraint, bool drawLeg=true) {

  const std::string leg1 = "w/o constraint";
//  const std::string leg2 = constraint == "topoConstr" ? "topo constraint" : constraint == "minvConstr" ? "m_{inv} constraint" : "ERROR";
  const std::string leg2 = "w/ constraint";

  TString currentMacroPath = __FILE__;
  TString directory = currentMacroPath(0, currentMacroPath.Last('/'));
  gROOT->Macro( directory + "/../styles/treeKF_qa2.dpg.style.cc" );

  TFile* fileIn1 = TFile::Open(fileName1.c_str());
  TFile* fileIn2 = TFile::Open(fileName2.c_str());
  if(fileIn1 == nullptr || fileIn2 == nullptr) throw std::runtime_error("fileIn1 == nullptr || fileIn2 == nullptr");

  TMultiGraph* mGr = new TMultiGraph();
  TGraphErrors* gr1 = fileIn1->Get<TGraphErrors>("sigma");
  TGraphErrors* gr2 = fileIn2->Get<TGraphErrors>("sigma");

  for(auto& g : {gr1, gr2}) {
    if(g == nullptr) throw std::runtime_error("g == nullptr");
    g->UseCurrentStyle();
    g->SetTitle("");
    g->SetMarkerStyle(kFullSquare);
    g->SetMarkerSize(3.2);
    g->SetLineWidth(3);
  }
  gr2->SetMarkerStyle(kOpenSquare); // TODO ad. hoc.
  gr1->SetMarkerColor(kBlue);
  gr2->SetMarkerColor(kRed);
  gr1->SetLineColor(kBlue);
  gr2->SetLineColor(kRed);
  SlightlyShiftXAxis(gr2);

  mGr->Add(gr1, "P");
  mGr->Add(gr2, "P");
  mGr->GetXaxis()->SetTitle(gr1->GetXaxis()->GetTitle());
  mGr->GetYaxis()->SetTitle(gr1->GetYaxis()->GetTitle());

  TLegend* leg = new TLegend(0.2, 0.7, 0.4, 0.9);
  leg->SetTextSize(0.07);
  leg->AddEntry(gr1, leg1.c_str(), "PE");
  leg->AddEntry(gr2, leg2.c_str(), "PE");

  TCanvas cc("cc", "cc", 1200, 800);
  mGr->Draw("AP");
  if(drawLeg) leg->Draw("same");
//  AddOneLineText(promptness, 0.74, 0.82, 0.87, 0.90);
  cc.Print("graph_qa3.pdf", "pdf");
}

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./graph_qa3 fileName1 fileName2 promptness constraint drawLeg=true" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileName1 = argv[1];
  const std::string fileName2 = argv[2];
  const std::string promptness = argv[3];
  const std::string constraint = argv[4];
  const bool drawLeg = argc > 5 ? string_to_bool(argv[5]) : true;

  graph_qa3(fileName1, fileName2, promptness, constraint, drawLeg);

  return 0;
}