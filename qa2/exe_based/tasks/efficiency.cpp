//
// Created by oleksii on 19.02.25.
//
#include "Helper.hpp"

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TString.h>

#include <iostream>

using namespace Helper;

void efficiency(const std::string& fileName, bool isSaveToRoot) {
  TString currentMacroPath = __FILE__;
  TString directory = currentMacroPath(0, currentMacroPath.Last('/'));
  gROOT->Macro( directory + "/../styles/mc_qa2.style.cc" );

  TFile* fileIn = TFile::Open(fileName.c_str());
  if(fileIn == nullptr) {
    throw std::runtime_error("fileIn == nullptr");
  }

  std::string fileOutName = "efficiency";

  struct Variable {
    std::string name_;
    int rebin_factor_;
  };

  std::vector<Variable> variables {
    {"Y",  2},
    {"Pt", 5},
    {"L",  5},
    {"T",  5},
  };

  TLegend* leg = new TLegend(0.2, 0.7, 0.4, 0.82);
  leg->SetBorderSize(0);

  TFile* fileOut{nullptr};
  if(isSaveToRoot) fileOut = TFile::Open((fileOutName + ".root").c_str(), "recreate");

  std::string printing_bracket = "(";
  for(auto& var : variables) {
    TH1D* histoGenPrompt = fileIn->Get<TH1D>(("EfficiencyQA/prompt/Generated_prompt/gen_" + var.name_).c_str());
    TH1D* histoRecPrompt = fileIn->Get<TH1D>(("EfficiencyQA/prompt/Candidates_prompt/rec_" + var.name_).c_str());
    TH1D* histoGenNonPrompt = fileIn->Get<TH1D>(("EfficiencyQA/nonprompt/Generated_nonprompt/gen_" + var.name_).c_str());
    TH1D* histoRecNonPrompt = fileIn->Get<TH1D>(("EfficiencyQA/nonprompt/Candidates_nonprompt/rec_" + var.name_).c_str());
    for(auto& histo : std::array<TH1D*, 4>{histoGenPrompt, histoRecPrompt, histoGenNonPrompt, histoRecNonPrompt}) {
      if(histo == nullptr) throw std::runtime_error("One of the TH1 histograms is nullptr!");
      if(var.rebin_factor_ != 1) histo->Rebin(var.rebin_factor_);
      histo->Sumw2();
      histo->UseCurrentStyle();
    }

    TH1D* histoEffPrompt = (TH1D*)histoRecPrompt->Clone();
    histoEffPrompt->Divide(histoGenPrompt);
    TH1D* histoEffNonPrompt = (TH1D*)histoRecNonPrompt->Clone();
    histoEffNonPrompt->Divide(histoGenNonPrompt);

    if(printing_bracket == "(") { // first canvas
      leg->AddEntry(histoEffPrompt, "prompt", "L");
      leg->AddEntry(histoEffNonPrompt, "nonprompt", "L");
    }

    for(auto& histo : std::array<TH1D*, 2>{histoEffPrompt, histoEffNonPrompt}) {
      histo->GetYaxis()->SetTitle("#varepsilon, %");
      histo->Scale(100);
      histo->SetTitle("");
    }
    histoEffPrompt->SetName(("hEffPrompt" + var.name_).c_str());
    histoEffNonPrompt->SetName(("hEffNonPrompt" + var.name_).c_str());
    histoEffPrompt->SetLineColor(kRed);
    histoEffNonPrompt->SetLineColor(kBlue);
    CustomizeHistogramsYRange({histoEffPrompt, histoEffNonPrompt}, 100);

    TCanvas ccEff("ccEff", "ccEff", 1200, 800);
    histoEffPrompt->Draw("");
    histoEffNonPrompt->Draw("same");
    leg->Draw("same");
    ccEff.Print((fileOutName + ".pdf" + printing_bracket).c_str(), "pdf");
    printing_bracket = "";
  } // variables

  TH2D* histoGenPrompt = fileIn->Get<TH2D>("EfficiencyQA/prompt/Generated_prompt/gen_Y_Vs_Pt");
  TH2D* histoRecPrompt = fileIn->Get<TH2D>("EfficiencyQA/prompt/Candidates_prompt/rec_Y_Vs_Pt");
  TH2D* histoGenNonPrompt = fileIn->Get<TH2D>("EfficiencyQA/nonprompt/Generated_nonprompt/gen_Y_Vs_Pt");
  TH2D* histoRecNonPrompt = fileIn->Get<TH2D>("EfficiencyQA/nonprompt/Candidates_nonprompt/rec_Y_Vs_Pt");
  for(auto& histo : std::array<TH2D*, 4>{histoGenPrompt, histoRecPrompt, histoGenNonPrompt, histoRecNonPrompt}) {
    if(histo == nullptr) throw std::runtime_error("One of the TH2 histograms is nullptr!");
    histo->Sumw2();
    histo->UseCurrentStyle();
  }

  TH2D* histoEffPrompt = (TH2D*)histoRecPrompt->Clone();
  histoEffPrompt->Divide(histoGenPrompt);
  TH2D* histoEffNonPrompt = (TH2D*)histoRecNonPrompt->Clone();
  histoEffNonPrompt->Divide(histoGenNonPrompt);

  histoEffPrompt->SetName("hEffPromptY_Vs_Pt");
  histoEffNonPrompt->SetName("hEffNonPromptY_Vs_Pt");

  for(auto& histo : std::array<TH2D*, 2>{histoEffPrompt, histoEffNonPrompt}) {
    histo->GetZaxis()->SetTitle("#varepsilon, %");
    histo->Scale(100);
    histo->SetMaximum(100);
    histo->SetTitle("");
    if(isSaveToRoot) {
      fileOut->cd();
      histo->Write();
    }
  }

  TCanvas ccEff2Prompt("ccEff2Prompt", "ccEff2Prompt", 1200, 800);
  ccEff2Prompt.SetRightMargin(0.16);
  histoEffPrompt->Draw("colz");
  AddOneLineText("prompt", 0.74, 0.82, 0.87, 0.90);
  ccEff2Prompt.Print((fileOutName + ".pdf").c_str(), "pdf");

  TCanvas ccEff2NonPrompt("ccEff2NonPrompt", "ccEff2NonPrompt", 1200, 800);
  ccEff2NonPrompt.SetRightMargin(0.16);
  histoEffNonPrompt->Draw("colz");
  AddOneLineText("nonprompt", 0.74, 0.82, 0.87, 0.90);
  ccEff2NonPrompt.Print((fileOutName + ".pdf)").c_str(), "pdf");

  if(isSaveToRoot) {
    fileOut->cd();
    histoEffPrompt->Write();
    histoEffNonPrompt->Write();
    fileOut->Close();
  }

  fileIn->Close();
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./efficiency fileName (isSaveRoot=false)" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileName = argv[1];
  const bool isSaveToRoot = argc > 2 ? string_to_bool(argv[2]) : false;

  efficiency(fileName, isSaveToRoot);

  return 0;
}