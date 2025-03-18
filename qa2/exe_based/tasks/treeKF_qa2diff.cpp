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

void treeKF_qa2diff(const std::string& fileName, int prompt_or_nonprompt, bool isDoFit, bool isSaveRoot) {
  TString currentMacroPath = __FILE__;
  TString directory = currentMacroPath(0, currentMacroPath.Last('/'));
  gROOT->Macro( directory + "/../styles/treeKF_qa2.dpg.style.cc" );

  TFile* fileIn = TFile::Open(fileName.c_str());
  if(fileIn == nullptr) {
    throw std::runtime_error("fileIn == nullptr");
  }

  std::string promptness;
  switch (prompt_or_nonprompt) {
    case 1: promptness = "prompt"; break;
    case 2: promptness = "nonprompt"; break;
    case -999: promptness = "data"; break;
    default: throw std::runtime_error("Illegal value of prompt_or_nonprompt. Should be either 1, 2 or -999");
  }

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

    auto cutsVar = FindCuts(fileIn, "TopoQA/Candidates_" + promptness + "_" + statusDcaFSel + "_" + var.cut_name_);

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
      const std::string histoName = "TopoQA/Candidates_"  + promptness + "_" + statusDcaFSel + "_" + cutName + "/Mass_"  + promptness + "_" + statusDcaFSel + "_" + cutName;
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
      TPaveText quant_text = ConvertHistoQuantitiesToText(quant, 0.70, 0.6, 0.90, 0.8);
      quant_text.Draw("same");

      if(isDoFit) {
        ShapeFitter shFtr(hIn);
        shFtr.SetExpectedMu(massLambdaC);
        shFtr.SetExpectedSigma(quant.stddev_);
        shFtr.SetPeakShape("DSCB");
        shFtr.SetBgPolN(-1); // FIXME effectively should zero BG polynomial. Careful, not tested yet! To be re-implemented properly.
        shFtr.Fit();
        TPaveText* fit_text = shFtr.ConvertFitParametersToText("peak", {0.20, 0.90});
        fit_text->Draw("same");
      }

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

    auto trueMassLine = HorizontalLine4Graph(massLambdaC, grMu);
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
    std::cout << " ./treeKF_qa2diff fileName (prompt_or_nonprompt=1 isDoFit=false isSaveRoot=false)" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileName = argv[1];
  const int prompt_or_nonprompt = argc>2 ? atoi(argv[2]) : 1;
  const bool isDoFit = argc > 3 ? string_to_bool(argv[3]) : false;
  const bool isSaveRoot = argc > 4 ? string_to_bool(argv[4]) : false;

  treeKF_qa2diff(fileName, prompt_or_nonprompt, isDoFit, isSaveRoot);

  return 0;
}


