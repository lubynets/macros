//
// Created by oleksii on 18.03.25.
//

#include "Helper.hpp"
#include "ShapeFitter.hpp"

#include <TCanvas.h>
#include <TROOT.h>

#include <iostream>

using namespace Helper;

void mass_fit(const std::string& fileName, bool isMC, bool isSaveToRoot) {
  TString currentMacroPath = __FILE__;
  TString directory = currentMacroPath(0, currentMacroPath.Last('/'));
  gROOT->Macro( directory + "/../styles/mc_qa2.style.cc" );

  TFile* fileIn = OpenFileWithNullptrCheck(fileName);

//  const std::vector<std::string> mcDataTypes {
//    "prompt",
//    "nonprompt",
//    "signal",
//    "background",
//  };

  const int rebinFactor = 4;

  const std::string mainDataType = isMC ? "all" : "data";

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

  for(auto& wu : weightsUsages) {
    std::string printingBracket = "(";
    for(auto& sc : sliceCuts) {
      const std::string histoName = mainDataType + "/Candidates_" + mainDataType + "_T_" + sc.first + "_" + sc.second + wu.dir_name_suffix_ + "/Mass" + wu.histo_name_suffix_ + "_all_T_" + sc.first + "_" + sc.second;
      const std::string cutRangeText = sc.first + " < T < " + sc.second + " (fs)";
      TH1D* histoIn = GetObjectWithNullptrCheck<TH1D>(fileIn, histoName);
      histoIn->UseCurrentStyle();
      histoIn->Sumw2();
      if(rebinFactor != 1) histoIn->Rebin(rebinFactor);

      ShapeFitter shapeFitter(histoIn);
      shapeFitter.SetExpectedMu(massLambdaC);
      shapeFitter.SetExpectedSigma(massLambdaCDetectorWidth);
      shapeFitter.SetSideBands(2.12, 2.20, 2.38, 2.42);
      shapeFitter.SetPeakShape("DSCB");
      shapeFitter.SetBgPolN(2);
      shapeFitter.Fit();

      TCanvas ccSideBand("ccSideBand", "");
      ccSideBand.SetCanvasSize(1200, 800);
      shapeFitter.GetSideBandHisto()->Draw();
      shapeFitter.GetSideBandFunc()->Draw("same");
      AddOneLineText(cutRangeText, {0.74, 0.82, 0.87, 0.90});
      TPaveText* parsSideBandText = shapeFitter.ConvertFitParametersToText("bg", {0.2, 0.9});
      parsSideBandText->Draw("same");
      ccSideBand.Print(("sideBand_" + wu.name_ + printingBracket).c_str(), "pdf");

      TCanvas ccPeak("ccPeak", "");
      ccPeak.SetCanvasSize(1200, 800);
      shapeFitter.GetPeakHisto()->Draw();
      shapeFitter.GetPeakFunc()->Draw("same");
      AddOneLineText(cutRangeText, {0.74, 0.82, 0.87, 0.90});
      TPaveText* parsPeakText = shapeFitter.ConvertFitParametersToText("peak", {0.2, 0.9});
      parsPeakText->Draw("same");
      ccPeak.Print(("peak_" + wu.name_ + printingBracket).c_str(), "pdf");

      TCanvas ccAll("ccAll", "");
      ccAll.SetCanvasSize(1200, 800);
      shapeFitter.GetAllHisto()->Draw();
      shapeFitter.GetAllFunc()->Draw("same");
      AddOneLineText(cutRangeText, {0.74, 0.82, 0.87, 0.90});
      TPaveText* parsAllText = shapeFitter.ConvertFitParametersToText("all", {0.2, 0.9});
      parsAllText->Draw("same");
      ccAll.Print(("all_" + wu.name_ + printingBracket).c_str(), "pdf");

      printingBracket = "";
    } // sliceCuts
    TCanvas emptyCanvas("emptyCanvas", "");
    emptyCanvas.SetCanvasSize(1200, 800);
    printingBracket = "]";
    for(auto& ccType : {"sideBand", "peak", "all"}) {
      emptyCanvas.Print((static_cast<std::string>(ccType) + "_" + wu.name_ + printingBracket).c_str(), "pdf");
    }
  } // weightsUsages
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