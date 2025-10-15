//
// Created by oleksii on 15.10.2025.
//
#include "HelperGeneral.hpp"
#include "HelperMath.hpp"
#include "HelperPlot.hpp"
#include "ShapeFitter.hpp"

#include <TCanvas.h>
#include <TLegend.h>

#include <iostream>

using namespace HelperGeneral;
using namespace HelperMath;
using namespace HelperPlot;

const std::vector<double> lifetimeRanges{0.2, 0.35, 0.5, 0.7, 0.9, 1.6};
const bool drawRefit{true};

const int rebinFactor{4};
const std::string peakShape{"DSCB"};
const int bgShape{2};

void mass_fit2(const std::string& fileName, const bool isSaveToRoot) {
  const std::string fitShape = peakShape + "pol" + std::to_string(bgShape);
  LoadMacro("styles/mc_qa2.style.cc");
  TFile* fileIn = OpenFileWithNullptrCheck(fileName);

  for(size_t iT=0, nTs=lifetimeRanges.size()-1; iT<nTs; ++iT) {
    const std::string priBra = EvaluatePrintingBracket(nTs, iT);
    const std::string histoInName = "data/pT_0_20/T_" + to_string_with_precision(lifetimeRanges.at(iT), 2) + "_" +  to_string_with_precision(lifetimeRanges.at(iT+1), 2) + "/hM_NPgt0.01";
    const std::string cutRangeText = to_string_with_precision(lifetimeRanges.at(iT), 2) + " < T < " + to_string_with_precision(lifetimeRanges.at(iT+1), 2) + " (ps)";
    TH1* histoIn = GetObjectWithNullptrCheck<TH1>(fileIn, histoInName);
    histoIn->UseCurrentStyle();
    histoIn->SetMarkerSize(0);
    histoIn->SetLineColor(kBlue);
    if(rebinFactor != 1) histoIn->Rebin(rebinFactor);

    ShapeFitter shapeFitter(histoIn);
    shapeFitter.SetExpectedMu(massLambdaC);
    shapeFitter.SetExpectedSigma(massLambdaCDetectorWidth);
    shapeFitter.SetSideBands(2.12, 2.23, 2.34, 2.42);
    shapeFitter.SetPeakShape(peakShape);
    shapeFitter.SetBgPolN(bgShape);
    shapeFitter.Fit();

    for(const auto& lines : {shapeFitter.GetAllFunc(), shapeFitter.GetAllReFunc(), shapeFitter.GetSideBandFunc(), shapeFitter.GetSideBandReFunc()}) lines->SetLineWidth(3);
    for(const auto& alls : {shapeFitter.GetAllFunc(), shapeFitter.GetAllReFunc()}) alls->SetLineColor(kRed);
    for(const auto& bgs : {shapeFitter.GetSideBandFunc(), shapeFitter.GetSideBandReFunc()}) bgs->SetLineColor(kGreen+2);
    for(const auto& prefits : {shapeFitter.GetAllFunc(), shapeFitter.GetSideBandFunc(), shapeFitter.GetPeakFunc()}) if(drawRefit) prefits->SetLineStyle(7);
    shapeFitter.GetPeakFunc()->SetFillColorAlpha(kBlue, 0.2);
    shapeFitter.GetPeakReFunc()->SetFillColorAlpha(kOrange+2, 0.2);
    for(const auto& sigs : {shapeFitter.GetPeakFunc(), shapeFitter.GetPeakReFunc()}) {
      sigs->SetLineWidth(0);
      sigs->SetFillStyle(1000);
    }

    TCanvas ccSideBand("ccSideBand", "");
    ccSideBand.SetCanvasSize(1200, 800);
    const double legBgUpperBorder = 0.9-(shapeFitter.GetSideBandFunc()->GetNpar()+1)*0.04;
    TLegend legBg(0.2, legBgUpperBorder - 0.1, 0.35, legBgUpperBorder);
    shapeFitter.GetSideBandHisto()->Draw("PE");
    legBg.AddEntry(shapeFitter.GetSideBandHisto(), "Sideband", "PE");
    legBg.AddEntry(shapeFitter.GetSideBandFunc(), "Sideband pre-fit", "L");
    shapeFitter.GetSideBandFunc()->Draw("same");
    if(drawRefit) {
      shapeFitter.GetSideBandReFunc()->Draw("same");
      legBg.AddEntry(shapeFitter.GetSideBandReFunc(), "Sideband re-fit", "L");
    }
    AddOneLineText(cutRangeText, {0.74, 0.82, 0.87, 0.90});
    TPaveText* parsSideBandText = shapeFitter.ConvertFitParametersToText("bg", {0.2, 0.9});
    parsSideBandText->Draw("same");
    legBg.Draw("same");
    ccSideBand.Print(("sideBand_" + fitShape + ".pdf" + priBra).c_str(), "pdf");

    TCanvas ccPeak("ccPeak", "");
    ccPeak.SetCanvasSize(1200, 800);
    const float legSigUpperBorder = 0.9-(shapeFitter.GetPeakFunc()->GetNpar()+1)*0.04;
    TLegend legSig(0.2, legSigUpperBorder - 0.1, 0.35, legSigUpperBorder);
    shapeFitter.GetPeakHisto()->SetMinimum(0);
    shapeFitter.GetPeakHisto()->Draw("PE");
    legSig.AddEntry(shapeFitter.GetPeakHisto(), "All - sideband fit", "PE");
    legSig.AddEntry(shapeFitter.GetPeakFunc(), "(All - sideband fit) fit", "F");
    shapeFitter.GetPeakFunc()->Draw("same");
    if(drawRefit) shapeFitter.GetPeakReFunc()->Draw("same");
    AddOneLineText(cutRangeText, {0.74, 0.82, 0.87, 0.90});
    TPaveText* parsPeakText = shapeFitter.ConvertFitParametersToText("peak", {0.2, 0.9});
    parsPeakText->Draw("same");
    legSig.Draw("same");
    ccPeak.Print(("peak_" + fitShape + ".pdf" + priBra).c_str(), "pdf");

    TCanvas ccAll("ccAll", "");
    ccAll.SetCanvasSize(1200, 800);
    shapeFitter.GetAllHisto()->SetMinimum(0);
    shapeFitter.GetAllHisto()->Draw("PE");
    shapeFitter.GetSideBandFunc()->Draw("same");
    if(drawRefit) shapeFitter.GetSideBandReFunc()->Draw("same");
    shapeFitter.GetAllFunc()->Draw("same");
    if(drawRefit) shapeFitter.GetAllReFunc()->Draw("same");
    shapeFitter.GetPeakFunc()->Draw("same FC");
    if(drawRefit) shapeFitter.GetPeakReFunc()->Draw("same FC");
    AddOneLineText(cutRangeText, {0.74, 0.82, 0.87, 0.90});
    TPaveText* parsAllText = shapeFitter.ConvertFitParametersToText("all", {0.2, 0.9});
    parsAllText->Draw("same");
    const float legAllUpperBorder = 0.9-(shapeFitter.GetAllFunc()->GetNpar()+1)*0.04;
    TLegend legAll(0.2, legAllUpperBorder - 0.2, 0.35, legAllUpperBorder);
    legAll.AddEntry(shapeFitter.GetAllFunc(), "Sig + Bg (pre-fit)", "L");
    legAll.AddEntry(shapeFitter.GetSideBandFunc(), "Bg (pre-fit)", "L");
    legAll.AddEntry(shapeFitter.GetPeakFunc(), "Sig (pre-fit)", "F");
    if(drawRefit) {
      legAll.AddEntry(shapeFitter.GetAllReFunc(), "Sig + Bg (re-fit)", "L");
      legAll.AddEntry(shapeFitter.GetSideBandReFunc(), "Bg (re-fit)", "L");
      legAll.AddEntry(shapeFitter.GetPeakReFunc(), "Sig (re-fit)", "F");
    }
    legAll.Draw("same");
    ccAll.Print(("all_" + fitShape + ".pdf" + priBra).c_str(), "pdf");
  } // nTs

  fileIn->Close();
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./mass_fit fileName (isSaveRoot=false)" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileName = argv[1];
  const bool isSaveToRoot = argc > 2 ? string_to_bool(argv[2]) : false;

  mass_fit2(fileName, isSaveToRoot);

  return 0;
}