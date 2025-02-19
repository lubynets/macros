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
    bool logy_eff_;
    bool logy_yield_;
  };

  std::vector<Variable> variables {
    {"Y",  2, false, false},
    {"Pt", 5, false, true },
    {"L",  5, false, true },
    {"T",  5, false, true },
  };

  TLegend* legEff = new TLegend(0.75, 0.7, 0.95, 0.82);
  TLegend* legYield = new TLegend(0.75, 0.7, 0.95, 0.82);
  legEff->SetBorderSize(0);
  legYield->SetBorderSize(0);

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
      legEff->AddEntry(histoEffPrompt, "prompt", "L");
      legEff->AddEntry(histoEffNonPrompt, "nonprompt", "L");
      legYield->AddEntry(histoRecPrompt, "rec", "L");
      legYield->AddEntry(histoGenPrompt, "gen", "L");
    }

    for(auto& histo : std::array<TH1D*, 2>{histoEffPrompt, histoEffNonPrompt}) {
      histo->GetYaxis()->SetTitle("#varepsilon, %");
      histo->Scale(100);
      histo->SetTitle("");
    }
    histoEffPrompt->SetName(("hEffPrompt" + var.name_).c_str());
    histoEffNonPrompt->SetName(("hEffNonPrompt" + var.name_).c_str());
    for(auto& histo : std::array<TH1D*, 3>{histoRecPrompt, histoRecNonPrompt, histoEffPrompt}) {
      histo->SetLineColor(kRed);
    }
    for(auto& histo : std::array<TH1D*, 3>{histoGenPrompt, histoGenNonPrompt, histoEffNonPrompt}) {
      histo->SetLineColor(kBlue);
    }
    CustomizeHistogramsYRange({histoRecPrompt, histoGenPrompt});
    CustomizeHistogramsYRange({histoRecNonPrompt, histoGenNonPrompt});
    CustomizeHistogramsYRange({histoEffPrompt, histoEffNonPrompt},100);

    auto PrintCanvas1D = [&] (const std::array<TH1D*, 2>& histos,
                            const std::string& namePrefix,
                            bool logy,
                            const std::string& oneLineText = "",
                            TLegend* leg = nullptr) {
        TCanvas cc("cc", "cc", 1200, 800);
        cc.SetLogy(logy);
        histos.at(0)->Draw("");
        histos.at(1)->Draw("same");
        if(logy) {
          for(auto& histo : histos) {
            histo->GetYaxis()->SetRangeUser(0.1, histo->GetMaximum()*2);
          }
        }
        if(leg != nullptr) leg->Draw("same");
        if(!oneLineText.empty()) AddOneLineText(oneLineText, 0.74, 0.82, 0.87, 0.90);
        cc.Print((namePrefix + ".pdf" + printing_bracket).c_str(), "pdf");
    };

    PrintCanvas1D({histoEffPrompt, histoEffNonPrompt}, "efficiency", var.logy_eff_, "", legEff);
    PrintCanvas1D({histoRecPrompt, histoGenPrompt}, "yieldPrompt", var.logy_yield_, "prompt", legYield);
    PrintCanvas1D({histoRecNonPrompt, histoGenNonPrompt}, "yieldNonPrompt", var.logy_yield_, "nonprompt", legYield);
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

  auto PrintCanvas2D = [&] (TH2D* histo,
                          const std::string& namePrefix,
                          const std::string& oneLineText="") {
      TCanvas cc("cc", "cc", 1200, 800);
      cc.SetRightMargin(0.16);
      histo->Draw("colz");
      if(!oneLineText.empty()) AddOneLineText(oneLineText, 0.69, 0.82, 0.82, 0.90);
      cc.Print((namePrefix + ".pdf" + printing_bracket).c_str(), "pdf");
  };

  PrintCanvas2D(histoEffPrompt, "efficiency", "prompt");
  PrintCanvas2D(histoRecPrompt, "yieldPrompt", "prompt, rec");
  PrintCanvas2D(histoRecNonPrompt, "yieldNonPrompt", "nonprompt, rec");
  printing_bracket = ")";
  PrintCanvas2D(histoEffNonPrompt, "efficiency", "nonprompt");
  PrintCanvas2D(histoGenPrompt, "yieldPrompt", "prompt, gen");
  PrintCanvas2D(histoGenNonPrompt, "yieldNonPrompt", "nonprompt, gen");

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