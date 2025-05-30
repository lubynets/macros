#include "Helper.hpp"
#include "ShapeFitter.hpp"

#include <TArrow.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TGraphMultiErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TPaveText.h>
#include <TROOT.h>

#include <vector>

using namespace Helper;

void mc_qa2diff(const std::string& fileName, int prompt_or_nonprompt, bool isDoFit, bool isSaveRoot) {
  TString currentMacroPath = __FILE__;
  TString directory = currentMacroPath(0, currentMacroPath.Last('/'));
  gROOT->Macro( directory + "/../styles/mc_qa2.dpg.style.cc" );

  if(prompt_or_nonprompt !=1 && prompt_or_nonprompt != 2) {
    throw std::runtime_error("prompt_or_nonprompt must be 1 or 2");
  }

  const std::string promptness = prompt_or_nonprompt == 1 ? "prompt" : "nonprompt";

  TFile* fileIn = TFile::Open(fileName.c_str());
  if(fileIn == nullptr) {
    throw std::runtime_error("fileIn == nullptr");
  }

  const std::string fileOutName = "mc_qa2diff";

  struct Variable {
    std::string name_;
    std::string cut_name_;
    std::string cut_title_;
    std::string cut_unit_;
    bool log_res_;
    bool log_pull_;
  };

  std::vector<Variable> vars {
//  name    cutname  title          unit    logres logpull
//    {"P",   "psim",  "#it{p}^{mc}",     "GeV/#it{c}", false, false},
    {"Pt",  "pTsim", "#it{p}_{T}^{mc}", "GeV/#it{c}", false, false},
//    {"Xpv", "nPVC",  "numPVTracks","",      false, false},
//    {"Ypv", "nPVC",  "numPVTracks","",      false, false},
//    {"Zpv", "nPVC",  "numPVTracks","",      false, false},
//    {"Xsv", "pTsim", "#it{p}_{T}^{mc}","GeV/#it{c}", false, false},
//    {"Ysv", "pTsim", "#it{p}_{T}^{mc}","GeV/#it{c}", false, false},
//    {"Zsv", "pTsim", "#it{p}_{T}^{mc}","GeV/#it{c}", false, false},
//    {"L",   "lsim",  "L^{mc}",     "cm",    false, false},
    {"T",   "tsim",  "T^{mc}",     "ps",    false, false},
//    {"T",   "trec",  "T^{rec}",    "ps",    false, false},
  };

  struct ResPull{
    std::string name_;
    std::string prefix_;
    bool is_draw_oneline_on_stddev_;
  };

  std::vector<ResPull> resPulls {
    {"Residual", "res",  false},
    {"Pull",     "pull", true}
  };

  for(auto& var : vars) {
    bool is_first_canvas{true};
    std::string printing_bracket = "(";
    const std::string xtitle = var.cut_unit_.empty() ? var.cut_title_ : var.cut_title_ + " (" + var.cut_unit_ + ")";

    auto cutsVar = FindCuts(fileIn, "PullsAndResiduals/CandidatesQA/noSel/"  + promptness + "/Candidates_Simulated_"  + promptness + "_total_" + var.cut_name_);

    std::vector<TGraphMultiErrors*> grMu;
    std::vector<TGraphErrors*> grSigma;
    grMu.resize(resPulls.size());
    grSigma.resize(resPulls.size());
    for(int iRP=0; iRP<resPulls.size(); iRP++) {
      grMu.at(iRP) = new TGraphMultiErrors(cutsVar.size(), 2);
      grSigma.at(iRP) = new TGraphErrors(cutsVar.size());
      for(auto& g : (std::vector<TGraph*>){grMu.at(iRP), grSigma.at(iRP)}) {
        g->SetTitle("");
        g->GetXaxis()->SetTitle(xtitle.c_str());
        g->SetMarkerStyle(kFullSquare);
        g->SetMarkerSize(3.2);
        g->SetMarkerColor(kBlue);
        g->SetLineWidth(3);
        g->SetLineColor(kBlue);
        g->SetFillColorAlpha(kBlue, 0.3);
      }
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

      for(int iRP=0; iRP<resPulls.size(); iRP++) {
        const std::string currentFileOutName = fileOutName + "_" + var.name_ + "_" +  resPulls.at(iRP).prefix_;
        const std::string cutName = var.cut_name_ + "_"  + cV.first + "_"  + cV.second;
        const std::string histoName = "PullsAndResiduals/CandidatesQA/noSel/"  + promptness + "/Candidates_Simulated_"  + promptness + "_total_"  + cutName + "/" +  resPulls.at(iRP).prefix_ + "_"  + var.name_ + "_"  + cutName;
        TH1D* hIn = fileIn->Get<TH1D>(histoName.c_str());
        if(hIn == nullptr) {
          throw std::runtime_error("hIn == nullptr for " + histoName);
        }

        if(is_first_canvas && iPoint==0) {
          const std::string residualSuffix = iRP == 0 ? " (" + var.cut_unit_ + ")" : "";
          grMu.at(iRP)->GetYaxis()->SetTitle((resPulls.at(iRP).name_ + " mean" + residualSuffix).c_str());
          grSigma.at(iRP)->GetYaxis()->SetTitle((resPulls.at(iRP).name_ + " width" + residualSuffix).c_str());
        }

        hIn->UseCurrentStyle();

        TCanvas cc(("cc" + resPulls.at(iRP).prefix_).c_str(), ("cc" + resPulls.at(iRP).prefix_).c_str(), 1200, 800);
        cc.SetLogy(var.log_res_);
        hIn->Draw();
        cutText.Draw("same");
        HistoQuantities quant = EvaluateHistoQuantities(hIn);
        TPaveText quant_text = ConvertHistoQuantitiesToText(quant, 0.70, 0.6, 0.90, 0.8);
        quant_text.Draw("same");

        if(isDoFit) {
          ShapeFitter shFtr(hIn);
          shFtr.SetExpectedMu(0);
          shFtr.SetExpectedSigma(quant.stddev_);
          shFtr.SetPeakShape("DoubleGaus");
          shFtr.SetBgPolN(-1); // FIXME effectively should zero BG polynomial. Careful, not tested yet! To be re-implemented properly.
          shFtr.Fit();
          TPaveText* fit_text = shFtr.ConvertFitParametersToText("peak", {0.20, 0.90});
          fit_text->Draw("same");
        }

        grMu.at(iRP)->SetPoint(iPoint, grX, quant.mean_);
        grMu.at(iRP)->SetPointEX(iPoint, grEX, grEX);
        grMu.at(iRP)->SetPointEY(iPoint, 0, quant.mean_err_, quant.mean_err_);
        grMu.at(iRP)->SetPointEY(iPoint, 1, quant.stddev_, quant.stddev_);

        grSigma.at(iRP)->SetPoint(iPoint, grX, quant.stddev_);
        grSigma.at(iRP)->SetPointError(iPoint, 0, quant.stddev_err_);

        if(is_first_canvas) printing_bracket = "(";
        else                printing_bracket = "";
        cc.Print((currentFileOutName + ".pdf" + printing_bracket).c_str(), "pdf");
      } // resPulls
      is_first_canvas = false;
      ++iPoint;
    } // cutsVar

    auto zeroLine = HorizontalLine4Graph(0, grMu.at(0));
    auto oneLine = HorizontalLine4Graph(1, grMu.at(0));
    SetLineDrawParameters({zeroLine, oneLine});

    TCanvas emptycanvas("emptycanvas", "", 1200, 800);
    for(int iRP=0; iRP<resPulls.size(); iRP++) {
      const std::string currentFileOutName = fileOutName + "_" + var.name_ + "_" +  resPulls.at(iRP).prefix_;
      TFile* fileOut{nullptr};
      if(isSaveRoot) fileOut = TFile::Open((currentFileOutName + ".root").c_str(), "recreate");

      TLegend legStat(0.24, 0.83, 0.37, 0.90);
      TLegend legBoth(0.24, 0.71, 0.37, 0.90);
      legStat.AddEntry(grMu.at(0), "stat. error", "E");
      legBoth.AddEntry(grMu.at(0), "stat. error", "E");
      auto entry = legBoth.AddEntry("", "distrib. width", "F");
      entry->SetFillColorAlpha(grMu.at(0)->GetFillColor(), 0.3);
      entry->SetLineColor(kWhite);
      entry->SetFillStyle(1000);

      TCanvas ccMuStat(("cc" +  resPulls.at(iRP).prefix_ + "Mu").c_str(), ("cc" +  resPulls.at(iRP).prefix_ + "Mu").c_str(), 1200, 800);
      grMu.at(iRP)->GetYaxis()->SetRangeUser(-0.199, 0.3);
      grMu.at(iRP)->Draw("AP; 2");
      zeroLine->Draw("L same");
//      legBoth.Draw("same");
//      if(var.name_.find("pv") == std::string::npos) AddOneLineText(promptness, 0.74, 0.82, 0.87, 0.90);
      ccMuStat.Print((currentFileOutName + ".pdf").c_str(), "pdf");
      if(isSaveRoot) grMu.at(iRP)->Write("mu");

      TCanvas ccMuWidth(("cc" +  resPulls.at(iRP).prefix_ + "Mu").c_str(), ("cc" +  resPulls.at(iRP).prefix_ + "Mu").c_str(), 1200, 800);
      grMu.at(iRP)->Draw("AP; X");
      if(iRP == 1 && var.name_ == "T") {
        for(int iOutlyerPoint=0; iOutlyerPoint<=1; iOutlyerPoint++) {
          TArrow* arrow = new TArrow(grMu.at(1)->GetPointX(iOutlyerPoint), 0.25, grMu.at(1)->GetPointX(iOutlyerPoint), 0.295, 0.01, "|>");
          arrow->SetLineColor(kBlue);
          arrow->Draw();
          AddOneLineText(to_string_with_significant_figures(grMu.at(1)->GetPointY(iOutlyerPoint), 2),
                         {static_cast<float>(grMu.at(1)->GetPointX(iOutlyerPoint)), 0.23,
                          static_cast<float>(grMu.at(1)->GetPointX(iOutlyerPoint)), 0.25}, "", 0.04);
        }
      }
//      CustomizeGraphYRange(grMu.at(iRP));
      zeroLine->Draw("L same");
//      legStat.Draw("same");
//      if(var.name_.find("pv") == std::string::npos) AddOneLineText(promptness, 0.74, 0.82, 0.87, 0.90);
      ccMuWidth.Print((currentFileOutName + ".pdf").c_str(), "pdf");

      TCanvas ccSigma(("cc" +  resPulls.at(iRP).prefix_ + "Sigma").c_str(), ("cc" +  resPulls.at(iRP).prefix_ + "Sigma").c_str(), 1200, 800);

      if(var.name_ == "T") grSigma.at(1)->GetXaxis()->SetLimits(0.01, 2);
      if(var.name_ == "Pt") grSigma.at(1)->GetXaxis()->SetLimits(0, 15.5);

      grSigma.at(iRP)->GetYaxis()->SetRangeUser(0.65, 1.3);
      grSigma.at(iRP)->Draw("AP");
      if(resPulls.at(iRP).is_draw_oneline_on_stddev_) oneLine->Draw("L same");
//      legStat.Draw("same");
//      if(var.name_.find("pv") == std::string::npos) AddOneLineText(promptness, 0.74, 0.82, 0.87, 0.90);
      ccSigma.Print((currentFileOutName + ".pdf").c_str(), "pdf");
      if(isSaveRoot) grSigma.at(iRP)->Write("sigma");

      emptycanvas.Print((currentFileOutName + ".pdf]").c_str(), "pdf");
      if(isSaveRoot) fileOut->Close();
    } // resPulls
  } // vars

  fileIn->Close();
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./mc_qa2diff fileName (prompt_or_nonprompt=1 isDoFit=true isSaveRoot=true)" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileName = argv[1];
  const int prompt_or_nonprompt = argc>2 ? atoi(argv[2]) : 1;
  const bool isDoFit = argc > 3 ? string_to_bool(argv[3]) : true;
  const bool isSaveRoot = argc > 4 ? string_to_bool(argv[4]) : true;

  mc_qa2diff(fileName, prompt_or_nonprompt, isDoFit, isSaveRoot);

  return 0;
}