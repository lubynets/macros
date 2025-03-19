//
// Created by oleksii on 18.03.25.
//

#include "Helper.hpp"
#include "ShapeFitter.hpp"

#include <TCanvas.h>
#include <TLegend.h>
#include <TROOT.h>

#include <iostream>

using namespace Helper;

std::vector<double> EvaluateLifetimeBinRanges(const std::vector<std::pair<std::string, std::string>>& sliceCuts, bool doPrint=false);

void mass_fit(const std::string& fileName, bool isMC, bool isSaveToRoot) {
  TString currentMacroPath = __FILE__;
  TString directory = currentMacroPath(0, currentMacroPath.Last('/'));
  gROOT->Macro( directory + "/../styles/mc_qa2.style.cc" );

  TFile* fileIn = OpenFileWithNullptrCheck(fileName);

  const int rebinFactor = 4;
//  const std::string peakShape = "DSCB";
//  const std::string peakShape = "Gaus";
  const std::string peakShape = "DoubleGaus";

  const int bgShape = 2;

  const bool drawRefit{false};
//  const bool drawRefit{true};

  const std::string mainDataType = isMC ? "all" : "data";

  std::string fitShape;
  if(peakShape == "DSCB") fitShape = "dscb.";
  if(peakShape == "Gaus") fitShape = "gaus.";
  if(peakShape == "DoubleGaus") fitShape = "doublegaus.";

  fitShape = fitShape + "pol" + std::to_string(bgShape);

  struct WeightsUsage {
    std::string name_;
    std::string dir_name_suffix_;
    std::string histo_name_suffix_;
  };

  std::vector<WeightsUsage> weightsUsages {
    {"wo_weight", "", ""},
    {"with_weight", "_weight_recEffWeigth", "_W"}
  };

  auto sliceCuts = FindCuts(fileIn, mainDataType + "/Candidates_" + mainDataType + "_T", true);
  auto ctBinEdges = EvaluateLifetimeBinRanges(sliceCuts);

  TFile* fileOut{nullptr};
  if(isSaveToRoot) fileOut = TFile::Open(("yields." + fitShape + ".root").c_str(), "recreate");

  for(auto& wu : weightsUsages) {
    std::string printingBracket = "(";
    TH1D* histoYieldSignal = new TH1D(("histoYieldSignal_" + wu.name_).c_str(), "", ctBinEdges.size()-1, ctBinEdges.data());
    histoYieldSignal->GetXaxis()->SetTitle("T (ps)");
    histoYieldSignal->GetYaxis()->SetTitle("Entries");
    int iSc{1};
    for(auto& sc : sliceCuts) {
      const std::string histoName = mainDataType + "/Candidates_" + mainDataType + "_T_" + sc.first + "_" + sc.second + wu.dir_name_suffix_ + "/Mass" + wu.histo_name_suffix_ + "_" + mainDataType + "_T_" + sc.first + "_" + sc.second;
      const std::string cutRangeText = sc.first + " < T < " + sc.second + " (ps)";
      TH1D* histoIn = GetObjectWithNullptrCheck<TH1D>(fileIn, histoName);
      histoIn->UseCurrentStyle();
      histoIn->Sumw2();
      if(rebinFactor != 1) histoIn->Rebin(rebinFactor);

      ShapeFitter shapeFitter(histoIn);
      shapeFitter.SetExpectedMu(massLambdaC);
      shapeFitter.SetExpectedSigma(massLambdaCDetectorWidth);
      shapeFitter.SetSideBands(2.12, 2.20, 2.38, 2.42);
      shapeFitter.SetPeakShape(peakShape);
      shapeFitter.SetBgPolN(bgShape);
      shapeFitter.Fit();

      histoYieldSignal->SetBinContent(iSc, shapeFitter.GetSignalIntegral3Sigma());
      histoYieldSignal->SetBinError(iSc, shapeFitter.GetSignalErrIntegral3Sigma());

      for(auto& lines : {shapeFitter.GetAllFunc(), shapeFitter.GetAllReFunc(), shapeFitter.GetSideBandFunc(), shapeFitter.GetSideBandReFunc()}) {
        lines->SetLineWidth(3);
      }
      for(auto& alls : {shapeFitter.GetAllFunc(), shapeFitter.GetAllReFunc()}) {
        alls->SetLineColor(kRed);
      }
      for(auto& bgs : {shapeFitter.GetSideBandFunc(), shapeFitter.GetSideBandReFunc()}) {
        bgs->SetLineColor(kGreen+2);
      }
      for(auto& prefits : {shapeFitter.GetAllFunc(), shapeFitter.GetSideBandFunc(), shapeFitter.GetPeakFunc()}) {
        if(drawRefit) prefits->SetLineStyle(7);
      }
      shapeFitter.GetPeakFunc()->SetFillColorAlpha(kBlue, 0.2);
      shapeFitter.GetPeakFunc()->SetLineWidth(0);
      shapeFitter.GetPeakFunc()->SetFillStyle(1000);

      TCanvas ccSideBand("ccSideBand", "");
      ccSideBand.SetCanvasSize(1200, 800);
      shapeFitter.GetSideBandHisto()->Draw();
      shapeFitter.GetSideBandFunc()->Draw("same");
      if(drawRefit) shapeFitter.GetSideBandReFunc()->Draw("same");
      AddOneLineText(cutRangeText, {0.74, 0.82, 0.87, 0.90});
      TPaveText* parsSideBandText = shapeFitter.ConvertFitParametersToText("bg", {0.2, 0.9});
      parsSideBandText->Draw("same");
      ccSideBand.Print(("sideBand_" + fitShape + "_" + wu.name_ + ".pdf" + printingBracket).c_str(), "pdf");

      TCanvas ccPeak("ccPeak", "");
      ccPeak.SetCanvasSize(1200, 800);
      shapeFitter.GetPeakHisto()->Draw();
      shapeFitter.GetPeakFunc()->Draw("same");
      if(drawRefit) shapeFitter.GetPeakReFunc()->Draw("same");
      AddOneLineText(cutRangeText, {0.74, 0.82, 0.87, 0.90});
      TPaveText* parsPeakText = shapeFitter.ConvertFitParametersToText("peak", {0.2, 0.9});
      parsPeakText->Draw("same");
      ccPeak.Print(("peak_" + fitShape + "_" + wu.name_ + ".pdf" + printingBracket).c_str(), "pdf");

      TCanvas ccAll("ccAll", "");
      ccAll.SetCanvasSize(1200, 800);
      shapeFitter.GetAllHisto()->Draw();
      shapeFitter.GetSideBandFunc()->Draw("same");
      if(drawRefit) shapeFitter.GetSideBandReFunc()->Draw("same");
      shapeFitter.GetAllFunc()->Draw("same");
      if(drawRefit) shapeFitter.GetAllReFunc()->Draw("same");
      shapeFitter.GetPeakFunc()->Draw("same FC");
      AddOneLineText(cutRangeText, {0.74, 0.82, 0.87, 0.90});
      TPaveText* parsAllText = shapeFitter.ConvertFitParametersToText("all", {0.2, 0.9});
      parsAllText->Draw("same");
      const float legendUpperBorder = 0.9-(shapeFitter.GetAllFunc()->GetNpar()+1)*0.04;
      TLegend leg(0.2, legendUpperBorder-0.1, 0.35, legendUpperBorder);
      leg.AddEntry(shapeFitter.GetAllFunc(), "Sig + Bg", "L");
      leg.AddEntry(shapeFitter.GetSideBandFunc(), "Bg", "L");
      leg.AddEntry(shapeFitter.GetPeakFunc(), "Sig", "F");
      leg.Draw("same");
      ccAll.Print(("all_" + fitShape + "_" + wu.name_ + ".pdf" + printingBracket).c_str(), "pdf");

      printingBracket = "";
      ++iSc;
    } // sliceCuts

    TCanvas ccHistoYield("ccHistoYield", "");
    ccHistoYield.SetCanvasSize(1200, 800);
    histoYieldSignal->Draw();
    ccHistoYield.Print(("all_" + fitShape + "_" + wu.name_ + ".pdf" + printingBracket).c_str(), "pdf");

    std::vector<std::string> ccS;
    for(auto& ccType : {"sideBand", "peak", "all"}) {
      ccS.emplace_back((static_cast<std::string>(ccType) + "_" + fitShape + "_" + wu.name_).c_str());
    }
    CloseCanvasPrinting(ccS);
    if(isSaveToRoot) {
      fileOut->cd();
      histoYieldSignal->Write();
    }
  } // weightsUsages
  if(isSaveToRoot) fileOut->Close();
}

std::vector<double> EvaluateLifetimeBinRanges(const std::vector<std::pair<std::string, std::string>>& sliceCuts, bool doPrint) {
  std::vector<double> result;
  for(auto& sc : sliceCuts) {
    result.emplace_back(atof(sc.first.c_str()));
  }
  result.emplace_back(atof(sliceCuts.back().second.c_str()));

  if(doPrint) {
    std::cout << "Info: EvaluateLifetimeBinRanges()\n";
    for(auto& r : result) {
      std::cout << r << "\t";
    }
    std::cout << "\n";
  }

  return result;
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./mass_fit fileName (isMc=true isSaveRoot=false)" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileName = argv[1];
  const bool isMc = argc > 2 ? string_to_bool(argv[2]) : true;
  const bool isSaveToRoot = argc > 3 ? string_to_bool(argv[3]) : false;

  mass_fit(fileName, isMc, isSaveToRoot);

  return 0;
}