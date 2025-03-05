//
// Created by oleksii on 19.02.25.
//
#include "Helper.hpp"

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLegendEntry.h>
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

//   const std::string selection = "isSel";
 const std::string selection = "noSel";

  const std::string fileOutName = "efficiency." + selection;

  struct Variable {
    std::string name_;
    int rebin_factor_;
    std::vector<double> rebin_edges_;
    bool logy_eff_;
    bool logy_yield_;
    std::string cut_text_;
  };

  std::vector<Variable> variables {
    {"Y",  2, {}, false, false, ""},
    {"Pt", 5, {}, false, true,  "|y| < 0.8" },
    {"L",  5, {}, false, true,  "|y| < 0.8"  },
    {"T",  5, {}, false, true,  "|y| < 0.8"  },
  };

  for(int iB=0; iB<40; iB++) {
    variables.at(3).rebin_edges_.push_back(0.025 * iB);
  }
  for(int iB=0; iB<11; iB++) {
    variables.at(3).rebin_edges_.push_back(1 + iB*0.1);
  }

  const std::vector<Color_t> colors{kRed, kBlue};

  TFile* fileOut{nullptr};
  if(isSaveToRoot) fileOut = TFile::Open((fileOutName + ".root").c_str(), "recreate");
  std::string printing_bracket = "(";

  for(auto& var : variables) {
    if(var.rebin_factor_ != 1 && !var.rebin_edges_.empty()) {
      std::cout << "Warning: For var " << var.name_ << " both rebin_factor_ and rebin_edges_ are set. Only rebin_edges_ will be applied\n";
    }
    TH1D* histoGenPrompt = fileIn->Get<TH1D>(("EfficiencyQA/prompt/gen/Generated_prompt/gen_" + var.name_).c_str());
    TH1D* histoRecPrompt = fileIn->Get<TH1D>(("EfficiencyQA/prompt/rec/" + selection + "/Candidates_prompt/rec_" + var.name_).c_str());
    TH1D* histoGenNonPrompt = fileIn->Get<TH1D>(("EfficiencyQA/nonprompt/gen/Generated_nonprompt/gen_" + var.name_).c_str());
    TH1D* histoRecNonPrompt = fileIn->Get<TH1D>(("EfficiencyQA/nonprompt/rec/" + selection + "/Candidates_nonprompt/rec_" + var.name_).c_str());
    histoGenPrompt->SetName((static_cast<std::string>(histoGenPrompt->GetName()) + "_Prompt").c_str());
    histoRecPrompt->SetName((static_cast<std::string>(histoRecPrompt->GetName()) + "_Prompt").c_str());
    histoGenNonPrompt->SetName((static_cast<std::string>(histoGenNonPrompt->GetName()) + "_NonPrompt").c_str());
    histoRecNonPrompt->SetName((static_cast<std::string>(histoRecNonPrompt->GetName()) + "_NonPrompt").c_str());
    auto ProcessHistogram = [&] (TH1D*& histo) {
      if(histo == nullptr) throw std::runtime_error("One of the TH1 histograms is nullptr!");
      histo->Sumw2();
      if(var.rebin_factor_ != 1 && var.rebin_edges_.empty()) histo->Rebin(var.rebin_factor_);
      if(!var.rebin_edges_.empty()) {
        std::cout << "Info: var " << var.name_ << " is rebinned to\n{";
        for(auto& be : var.rebin_edges_) {
          std::cout << be << " ";
        }
        std::cout << "}\n";
        histo = dynamic_cast<TH1D*>(histo->Rebin(var.rebin_edges_.size() - 1,histo->GetName(),var.rebin_edges_.data()));
      }
      histo->UseCurrentStyle();
    };
    ProcessHistogram(histoGenPrompt);
    ProcessHistogram(histoRecPrompt);
    ProcessHistogram(histoGenNonPrompt);
    ProcessHistogram(histoRecNonPrompt);

    TH1D* histoEffPrompt = (TH1D*)histoRecPrompt->Clone();
    histoEffPrompt->Divide(histoGenPrompt);
    TH1D* histoEffNonPrompt = (TH1D*)histoRecNonPrompt->Clone();
    histoEffNonPrompt->Divide(histoGenNonPrompt);

    for(auto& histo : std::array<TH1D*, 2>{histoEffPrompt, histoEffNonPrompt}) {
      histo->GetYaxis()->SetTitle("#varepsilon, %");
      histo->Scale(100);
      histo->SetTitle("");
    }
    histoEffPrompt->SetName(("hEffPrompt" + var.name_).c_str());
    histoEffNonPrompt->SetName(("hEffNonPrompt" + var.name_).c_str());

    auto PrintCanvas1D = [&] (std::vector<TH1*> histos,
                              const std::vector<Color_t >& colors,
                              const std::string& namePrefix,
                              bool logy,
                              const std::string& oneLineText = "",
                              const std::vector<std::string>& legTexts={}) {
      TLegend* leg = new TLegend(0.25, 0.6, 0.45, 0.72);
      leg->SetBorderSize(0);
      const double yAxisMin = 0;
      const double yAxisMax = static_cast<TString>(histos.at(0)->GetName()).Contains("Eff") ? 100 : 1e9;
      if(histos.size() != colors.size()) throw std::runtime_error("PrintCanvas1D(): histos.size() != colors.size()");
      for(int iH=0; iH<histos.size(); iH++) {
        histos.at(iH)->SetLineColor(colors.at(iH));
        leg->AddEntry(histos.at(iH), legTexts.at(iH).c_str(), "L");
      }
      CustomizeHistogramsYRange(histos, yAxisMin, yAxisMax);
      TCanvas cc("cc", "cc", 1200, 800);
      cc.SetLogy(logy);
      histos.at(0)->Draw("");
      histos.at(1)->Draw("same");
      if(logy) {
        for(auto& histo : histos) {
          histo->GetYaxis()->SetRangeUser(std::max(0.1, histo->GetMinimum()/2), histo->GetMaximum()*2);
        }
      }
      leg->Draw("same");
      AddOneLineText(oneLineText, 0.74, 0.82, 0.87, 0.90);
      AddOneLineText(var.cut_text_, 0.60, 0.82, 0.73, 0.90);
      cc.Print((namePrefix + ".pdf" + printing_bracket).c_str(), "pdf");
    };

    if(isSaveToRoot) {
      fileOut->cd();
      for(auto& histo : std::vector<TH1D*>{histoEffPrompt, histoEffNonPrompt, histoRecPrompt, histoRecNonPrompt, histoGenPrompt, histoGenNonPrompt}) {
        histo->Write();
      }
    }

    PrintCanvas1D({histoEffPrompt, histoEffNonPrompt}, colors, "efficiency." + selection, var.logy_eff_, "", {"prompt", "nonprompt"});
    PrintCanvas1D({histoRecPrompt, histoGenPrompt}, colors, "yieldPrompt." + selection, var.logy_yield_, "prompt", {"rec", "gen"});
    PrintCanvas1D({histoRecNonPrompt, histoGenNonPrompt}, colors, "yieldNonPrompt." + selection, var.logy_yield_, "nonprompt", {"rec", "gen"});
    PrintCanvas1D({histoRecPrompt, histoRecNonPrompt}, colors, "yieldRec." + selection, var.logy_yield_, "rec", {"prompt", "nonprompt"});
    PrintCanvas1D({histoGenPrompt, histoGenNonPrompt}, colors, "yieldGen." + selection, var.logy_yield_, "gen", {"prompt", "nonprompt"});
    printing_bracket = "";
  } // variables

  TH2D* histoGenPrompt = fileIn->Get<TH2D>("EfficiencyQA/prompt/gen/Generated_prompt/gen_Y_Vs_Pt");
  TH2D* histoRecPrompt = fileIn->Get<TH2D>(("EfficiencyQA/prompt/rec/" + selection + "/Candidates_prompt/rec_Y_Vs_Pt").c_str());
  TH2D* histoGenNonPrompt = fileIn->Get<TH2D>("EfficiencyQA/nonprompt/gen/Generated_nonprompt/gen_Y_Vs_Pt");
  TH2D* histoRecNonPrompt = fileIn->Get<TH2D>(("EfficiencyQA/nonprompt/rec/" + selection + "/Candidates_nonprompt/rec_Y_Vs_Pt").c_str());
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

  PrintCanvas2D(histoEffPrompt, "efficiency." + selection, "prompt");
  PrintCanvas2D(histoRecPrompt, "yieldPrompt." + selection, "prompt, rec");
  PrintCanvas2D(histoRecNonPrompt, "yieldNonPrompt." + selection, "nonprompt, rec");
  PrintCanvas2D(histoRecPrompt, "yieldRec." + selection, "prompt, rec");
  PrintCanvas2D(histoGenPrompt, "yieldGen." + selection, "prompt, gen");
  printing_bracket = ")";
  PrintCanvas2D(histoEffNonPrompt, "efficiency." + selection, "nonprompt");
  PrintCanvas2D(histoGenPrompt, "yieldPrompt." + selection, "prompt, gen");
  PrintCanvas2D(histoGenNonPrompt, "yieldNonPrompt." + selection, "nonprompt, gen");
  PrintCanvas2D(histoRecNonPrompt, "yieldRec." + selection, "nonprompt, rec");
  PrintCanvas2D(histoGenNonPrompt, "yieldGen." + selection, "nonprompt, gen");

  if(isSaveToRoot) fileOut->Close();
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
