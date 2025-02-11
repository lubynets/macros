//
// Created by oleksii on 10.02.25.
//
#include "Helper.hpp"

#include <TCanvas.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TROOT.h>

#include <iostream>

using namespace Helper;

void graph_qa3(const std::string& fileName1, const std::string& fileName2) {

  const std::string leg1 = "no constraint";
  const std::string leg2 = "m_{inv} constraint";
//  const std::string leg2 = "topo constraint";

//  const std::string promptness = "prompt";
  const std::string promptness = "nonprompt";

  TString currentMacroPath = __FILE__;
  TString directory = currentMacroPath(0, currentMacroPath.Last('/'));
  gROOT->Macro( directory + "/../styles/treeKF_qa2.style.cc" );

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
    g->SetMarkerSize(2);
    g->SetLineWidth(3);
  }
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
  leg->AddEntry(gr1, leg1.c_str(), "PE");
  leg->AddEntry(gr2, leg2.c_str(), "PE");

  TCanvas cc("cc", "cc", 1200, 800);
  mGr->Draw("AP");
  leg->Draw("same");
  AddOneLineText(promptness, 0.74, 0.82, 0.87, 0.90);
  cc.Print("graph_qa3.pdf", "pdf");
}

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./graph_qa3 fileName1 fileName2" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileName1 = argv[1];
  const std::string fileName2 = argv[2];

  graph_qa3(fileName1, fileName2);

  return 0;
}