//
// Created by oleksii on 31.01.25.
//

#include "Helper.hpp"
#include "ShapeFitter.hpp"

#include <TCanvas.h>
#include <TGraphMultiErrors.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TROOT.h>

#include <iostream>

using namespace Helper;

void treeKF_qa2diff(const std::string& fileName, int prompt_or_nonprompt, bool isSaveRoot) {
  TString currentMacroPath = __FILE__;
  TString directory = currentMacroPath(0, currentMacroPath.Last('/'));
  gROOT->Macro( directory + "/treeKF_qa2.style.cc" );

  TFile* fileIn = TFile::Open(fileName.c_str());
  if(fileIn == nullptr) {
    throw std::runtime_error("fileIn == nullptr");
  }

  const std::string promptness = prompt_or_nonprompt == 1 ? "prompt" : "nonprompt";

  const std::string statusDcaFSel = "isDcaFSel";
//  const std::string statusDcaFSel = "noDcaFSel";

  const std::string fileOutName = "treeKF_qa2diff";

  struct Variable {
    std::string name_;
    std::string cut_name_;
    std::string cut_title_;
    std::string cut_unit_;
    bool logy_;
  };

  std::vector<Variable> vars {
    {"Mass", "pT", "p_{T}", "GeV/c", false},
  };

  for(auto& var : vars) {
    bool is_first_canvas{true};
    std::string printing_bracket = "(";
    const std::string currentFileOutName = fileOutName + "_" + var.name_;
    TFile* fileOut{nullptr};
    if(isSaveRoot) fileOut = TFile::Open((currentFileOutName + ".root").c_str(), "recreate");

    auto cutsVar = FindCuts(fileIn, "Candidates_" + promptness + "_" + statusDcaFSel + "_" + var.cut_name_);

    TGraphMultiErrors* grMu{nullptr};
    TGraphErrors* grSigma{nullptr};

    grMu = new TGraphMultiErrors(cutsVar.size(), 2);
    grSigma = new TGraphErrors(cutsVar.size());
    for(auto& g : (std::vector<TGraph*>){grMu, grSigma}) {
      g->SetTitle("");
      g->GetXaxis()->SetTitle((var.cut_title_ + ", " + var.cut_unit_).c_str());
      g->SetMarkerStyle(kFullSquare);
      g->SetMarkerSize(2);
      g->SetMarkerColor(kBlue);
      g->SetLineWidth(3);
      g->SetLineColor(kBlue);
      g->SetFillColorAlpha(kBlue, 0.3);
    }

    TPaveText generalText(0.74, 0.82, 0.87, 0.90, "brNDC");
    generalText.SetFillColor(0);
    generalText.SetTextSize(0.03);
    generalText.SetTextFont(62);
    generalText.AddText(promptness.c_str());

    int iPoint{0};
    for(auto& cV : cutsVar) {
      TPaveText cutText = generalText;
      cutText.AddText((cV.first + " < " + var.cut_title_ + " < " + cV.second + " " + var.cut_unit_).c_str());
      const float grX = (atof(cV.first.c_str()) + atof(cV.second.c_str())) / 2;
      const float grEX = -(atof(cutsVar.at(0).second.c_str()) - atof(cutsVar.at(0).first.c_str())) / 15;

      const std::string cutName = var.cut_name_ + "_"  + cV.first + "_"  + cV.second;
      const std::string histoName = "Candidates_"  + promptness + "_" + statusDcaFSel + "_" + cutName + "/Mass_"  + promptness + "_" + statusDcaFSel + "_" + cutName;
      TH1D* hIn = fileIn->Get<TH1D>(histoName.c_str());
      if(hIn == nullptr) {
        throw std::runtime_error("hIn == nullptr for " + histoName);
      }

      if(is_first_canvas && iPoint==0) {
        grMu->GetYaxis()->SetTitle((std::string("Mean of ") + hIn->GetXaxis()->GetTitle()).c_str());
        grSigma->GetYaxis()->SetTitle((std::string("Width of ") + hIn->GetXaxis()->GetTitle()).c_str());
      }

      hIn->UseCurrentStyle();

      TCanvas cc("cc", "cc", 1200, 800);
      cc.SetLogy(var.logy_);
      hIn->Draw();
      cutText.Draw("same");
      HistoQuantities quant = EvaluateHistoQuantities(hIn);
      TPaveText quant_text = ConvertHistoQuantitiesToText(quant, 0.74, 0.6, 0.87, 0.7);
      quant_text.Draw("same");

      ShapeFitter shFtr(hIn);
      shFtr.SetExpectedMu(massLambdaC);
      shFtr.SetExpectedSigma(quant.stddev_);
      shFtr.Fit("DoubleGaus");
      TPaveText fit_text = shFtr.FitParametersToText(0.30, 0.79, 0.45, 0.89);

      grMu->SetPoint(iPoint, grX, quant.mean_);
      grMu->SetPointEX(iPoint, grEX, grEX);
      grMu->SetPointEY(iPoint, 0, quant.mean_err_, quant.mean_err_);
      grMu->SetPointEY(iPoint, 1, quant.stddev_, quant.stddev_);

      grSigma->SetPoint(iPoint, grX, quant.stddev_);
      grSigma->SetPointError(iPoint, 0, quant.stddev_err_);

      if(is_first_canvas) printing_bracket = "(";
      else                printing_bracket = "";
      cc.Print((currentFileOutName + ".pdf" + printing_bracket).c_str(), "pdf");
      is_first_canvas = false;
      ++iPoint;
    } // cutsVar

    auto trueMassLine = HorizontalLine4Graph(2.28646, grMu);
    SetLineDrawParameters({trueMassLine});

    TCanvas emptycanvas("emptycanvas", "", 1200, 800);
    TLegend legStat(0.24, 0.83, 0.37, 0.90);
    TLegend legBoth(0.24, 0.71, 0.37, 0.90);
    legStat.AddEntry(grMu, "stat. error", "E");
    legBoth.AddEntry(grMu, "stat. error", "E");
    auto entry = legBoth.AddEntry("", "distrib. width", "F");
    entry->SetFillColorAlpha(grMu->GetFillColor(), 0.3);
    entry->SetLineColor(kWhite);
    entry->SetFillStyle(1000);

    TCanvas ccMuStat("ccMuStat", "ccMuStat", 1200, 800);
    grMu->Draw("AP; 2");
    CustomizeGraphYRange(grMu, 2, trueMassLine);
    trueMassLine->Draw("L same");
    legBoth.Draw("same");
    ccMuStat.Print((currentFileOutName +".pdf").c_str(), "pdf");
    if(isSaveRoot) grMu->Write("mu");

    TCanvas ccMuWidth("ccMuWidth", "ccMuWidth", 1200, 800);
    grMu->Draw("AP; X");
    CustomizeGraphYRange(grMu, 1, trueMassLine);
    trueMassLine->Draw("L same");
    legStat.Draw("same");
    ccMuWidth.Print((currentFileOutName +".pdf").c_str(), "pdf");

    TCanvas ccSigma("ccSigma", "ccSigma", 1200, 800);
    grSigma->Draw("AP");
    legStat.Draw("same");
    ccSigma.Print((currentFileOutName +".pdf").c_str(), "pdf");
    if(isSaveRoot) grSigma->Write("sigma");

    emptycanvas.Print((currentFileOutName +".pdf]").c_str(), "pdf");
    if(isSaveRoot) fileOut->Close();
  } // vars
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./treeKF_qa2diff fileName (prompt_or_nonprompt isSaveRoot)" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileName = argv[1];
  const int prompt_or_nonprompt = argc>2 ? atoi(argv[2]) : 1;
  const bool isSaveRoot = argc >= 4 && strcmp(argv[3], "true") == 0;

  treeKF_qa2diff(fileName, prompt_or_nonprompt, isSaveRoot);

  return 0;
}


