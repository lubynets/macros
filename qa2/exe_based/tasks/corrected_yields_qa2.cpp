//
// Created by oleksii on 26.06.25.
//
#include "HelperCtFit.hpp"
#include "HelperGeneral.hpp"
#include "HelperMath.hpp"
#include "HelperPlot.hpp"

#include <TGraphErrors.h>
#include <TFile.h>
#include <TLegend.h>
#include <TLine.h>
#include <TNamed.h>

#include <boost/program_options.hpp>

#include <cmath>
#include <iostream>
#include <string>

using namespace HelperGeneral;
using namespace HelperMath;
using namespace HelperPlot;

constexpr bool IsSaveCanvasAsRoot{true};

constexpr HelperCtFit::RunModes RunMode {
//   HelperCtFit::RunModes::MeanFitOnly
  HelperCtFit::RunModes::AllWoOne
//   HelperCtFit::RunModes::AllPossible
};

enum UncModes {
  StatOnly = 0,
  SystOnly,
  StatAndSyst
};
constexpr int UncMode{StatOnly};

void corrected_yields_qa2(const std::string& fileNameCutVar, const std::string& fileNameMC, const std::string& fileNameSystUnc) {
  LoadMacro("styles/mc_qa2.style.cc");

  const bool isMc = !fileNameMC.empty();
  const bool isSystUnc = UncMode > StatOnly;

  if(isSystUnc && fileNameSystUnc.empty()) throw std::runtime_error("corrected_yields_qa2() systematic uncertainties file is not provided");

  const std::vector<double> lifetimeRanges = {0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8};
  const std::string integralOption = "I";

  struct Promptness {
    std::string name_;
    std::string histo_name_;
  };

  const std::vector<Promptness> promptnesses {
    {"prompt",    "Prompt"},
    {"nonprompt", "NonPrompt"}
  };

  TFile* fileCutVar = OpenFileWithNullptrCheck(fileNameCutVar);
  TFile* fileMC = isMc ? OpenFileWithNullptrCheck(fileNameMC) : nullptr;
  TFile* fileSystUnc = isSystUnc ? OpenFileWithNullptrCheck(fileNameSystUnc) : nullptr;

  TH1* histoSystUnc = isSystUnc ? GetObjectWithNullptrCheck<TH1>(fileSystUnc, "systErrorsGet") : nullptr;
  if(isSystUnc) {
    CheckTAxisForRanges(*histoSystUnc->GetXaxis(), lifetimeRanges);
    histoSystUnc = dynamic_cast<TH1D*>(histoSystUnc->Rebin(lifetimeRanges.size() - 1, histoSystUnc->GetName(), lifetimeRanges.data()));
  }

  TH1* histoStatus{};
  histoStatus = fileCutVar->Get<TH1>("hMinimizationStatus");
  int bestDropSet{HelperGeneral::UndefValueInt};
  if (histoStatus != nullptr) {
    CheckTAxisForRanges(*histoStatus->GetXaxis(), lifetimeRanges);
    histoStatus = dynamic_cast<TH1D*>(histoStatus->Rebin(lifetimeRanges.size() - 1, histoStatus->GetName(), lifetimeRanges.data()));
    std::vector<int> badBins{};
    for(int iBin=1, nBins=histoStatus->GetNbinsX(); iBin<=nBins; ++iBin) {
      if(histoStatus->GetBinContent(iBin) != 1) badBins.push_back(iBin);
    }
    bestDropSet = HelperCtFit::EvalDropSet(badBins);
  }
  std::cout << "bestDropSet = " << bestDropSet << "\n";
  TNamed bestDropSetName("bestDropSet", std::to_string(bestDropSet));

  for (size_t iP = 0, nP = promptnesses.size(); iP < nP; ++iP) {
    const std::string promptness = promptnesses.at(iP).name_;

    TH1* hCutVar = GetObjectWithNullptrCheck<TH1>(fileCutVar, "hCorrYields" + promptnesses.at(iP).histo_name_);
    CheckTAxisForRanges(*hCutVar->GetXaxis(), lifetimeRanges);
    hCutVar = dynamic_cast<TH1D*>(hCutVar->Rebin(lifetimeRanges.size() - 1, hCutVar->GetName(), lifetimeRanges.data()));
    hCutVar->UseCurrentStyle();
    hCutVar->SetLineColor(kRed);
    hCutVar->SetMarkerColor(kRed);

    if(promptness == "prompt" && isSystUnc) {
      for(int iCt=1, nCts=hCutVar->GetNbinsX(); iCt<=nCts; ++iCt) {
        const double stat = hCutVar->GetBinError(iCt);
        const double syst = histoSystUnc->GetBinContent(iCt);
        const double unc = UncMode == SystOnly ? syst : std::sqrt(stat*stat + syst*syst);
        hCutVar->SetBinError(iCt, unc);
      }
    }

    TH1* hMC = isMc ? GetObjectWithNullptrCheck<TH1>(fileMC, "gen/" + promptness + "/hT") : nullptr;
    if (isMc) {
      CheckTAxisForRanges(*hMC->GetXaxis(), lifetimeRanges);
      hMC->UseCurrentStyle();
      hMC->SetLineColor(kBlue);
      hMC->SetMarkerColor(kBlue);

      hMC = dynamic_cast<TH1D*>(hMC->Rebin(lifetimeRanges.size() - 1, hMC->GetName(),lifetimeRanges.data()));

      CustomizeHistogramsYRange({hCutVar, hMC}, true);
    } // isMc

    auto PrintYields = [&](const std::string& ccName, TH1* hCutVar, TH1* hMc, const std::string& promptness, const std::string& priBra) {
      TCanvas cc("cc", "");
      cc.SetLogy();
      cc.SetCanvasSize(1200, 800);
      if (isMc) hMc->Draw("E1");
      hCutVar->Draw("E1 same");
      const double legy1 = isMc ? 0.7 : 0.8;
      TLegend legPrompt(0.7, legy1, 0.9, 0.9);
      legPrompt.SetHeader(promptness.c_str());
      if (isMc) legPrompt.AddEntry(hMc, "MC", "PL");
      legPrompt.AddEntry(hCutVar, "Cut Var", "PL");
      legPrompt.Draw("same");
      cc.Print((ccName + ".pdf" + priBra).c_str(), "pdf");
    };
    if(isMc) PrintYields("compar", hCutVar, hMC, promptness, EvaluatePrintingBracket(promptnesses, iP));

    if(promptness != "prompt") continue;

    TH1* hCutVarDiff = dynamic_cast<TH1*>(hCutVar->Clone());
    TH1* hMcDiff = isMc ? dynamic_cast<TH1*>(hMC->Clone()) : nullptr;
    for(const auto& h : {hCutVarDiff, hMcDiff}) {
      if(h == nullptr) continue;
      h->Scale(1.0, "width");
      h->GetXaxis()->SetTitle("#it{t} (ps)");
      h->GetYaxis()->SetTitle("d#it{N}/d#it{t} (ps^{-1})");
    }

    auto graphsQa = HelperCtFit::CreateQaGraphs();

    const int nBins = hCutVarDiff->GetNbinsX();
    const int nDropSets = RunMode > HelperCtFit::RunModes::MeanFitOnly ? 1 << nBins : 1;
    double referenceSigmaTau{HelperGeneral::UndefValueDouble};

    TFile* fileOut = IsSaveCanvasAsRoot ? TFile::Open("ctfit.root", "recreate") : nullptr;
    if(IsSaveCanvasAsRoot) bestDropSetName.Write();

    TCanvas emptycanvas("emptycanvas", "", 1200, 800);
    emptycanvas.Print("ctfit.pdf[", "pdf");

    auto FitLifetimeHistos = [&](TH1* histoRec, TH1* histoMc, int dropSet) {
      const auto binsToDrop = HelperCtFit::EvalBinsToDrop(dropSet);

      TF1* fitCutVar = FitLifetimeHisto(histoRec, integralOption);
      TF1* fitMc = isMc ? FitLifetimeHisto(histoMc, integralOption) : nullptr;

      if(dropSet == 0) referenceSigmaTau = fitCutVar->GetParError(1);
      HelperCtFit::FillQaGraphs(graphsQa, fitCutVar, dropSet, referenceSigmaTau);

      const std::pair<std::string, std::string> fitResultsCutVar = HelperCtFit::FitResults(fitCutVar);
      const std::pair<std::string, std::string> fitResultsMC = isMc ? HelperCtFit::FitResults(fitMc, "[MC]") : std::pair<std::string, std::string>();

      const float textX1 = 0.70;
      const float textX2 = 0.86;
      const float textY2 = 0.86;
      const float textYStep = 0.06;

      TCanvas ccFit("ccFit", "");
      ccFit.SetCanvasSize(1200, 800);
      ccFit.SetLogy();
      if (isMc) {
        histoMc->Draw("E1");
        fitMc->Draw("same");
        AddOneLineText(fitResultsMC.first, {textX1, textY2 - 1*textYStep, textX2, textY2}, "brNDC", 0.04);
        AddOneLineText(fitResultsMC.second, {textX1, textY2 - 2*textYStep, textX2, textY2 - 1*textYStep}, "brNDC", 0.04);
      }
      histoRec->Draw("E1 same");

      std::string droppedBins{};
      for(const auto& binToDrop : binsToDrop) {
        droppedBins.append(std::to_string(binToDrop));
        droppedBins.append(", ");
      }
      if(dropSet != 0) {
        droppedBins.pop_back();
        droppedBins.pop_back();
        histoRec->SetTitle(("drop set = " + std::to_string(dropSet) + ", dropped bins: {" + droppedBins + "}").c_str());
      }
      fitCutVar->Draw("same");
      AddOneLineText(fitResultsCutVar.first, {textX1, textY2 - 3*textYStep, textX2, textY2 - 2*textYStep}, "brNDC", 0.04);
      AddOneLineText(fitResultsCutVar.second, {textX1, textY2 - 4*textYStep, textX2, textY2 - 3*textYStep}, "brNDC", 0.04);
      AddOneLineText(HelperCtFit::LifetimePdgSv, {textX1, textY2 - 6*textYStep, textX2, textY2 - 5*textYStep}, "brNDC", 0.04);
      ccFit.Print("ctfit.pdf", "pdf");

      if(IsSaveCanvasAsRoot && dropSet == 0) {
        fileOut->cd();
        ccFit.Write();
      }
      const std::string fitCutVarSaveName = dropSet == 0 ? "fitFunc" : "fitFunc_dropSet" + std::to_string(dropSet);
      if(IsSaveCanvasAsRoot) fitCutVar->Write(fitCutVarSaveName.c_str());

      TH1* hCutVarRatio = dynamic_cast<TH1*>(histoRec->Clone());
      TH1* hMcRatio = isMc ? dynamic_cast<TH1*>(histoMc->Clone()) : nullptr;
      TH1* hCutVarResidual = dynamic_cast<TH1*>(histoRec->Clone());
      TH1* hMcResidual = isMc ? dynamic_cast<TH1*>(histoMc->Clone()) : nullptr;
      ScalePlotVertically(hCutVarRatio, histoRec, 2);
      hCutVarRatio->GetYaxis()->SetTitle("Data / Fit");
      DivideHistoByFunction(hCutVarRatio, fitCutVar, integralOption);
      ScalePlotVertically(hCutVarResidual, histoRec, 2);
      hCutVarResidual->GetYaxis()->SetTitle("(Data - Fit) / #sigma_{Data}");
      EvalNormDifferenceHistoFromFunction(hCutVarResidual, fitCutVar, integralOption);
      CustomizeHistogramsYRange({hCutVarRatio}, false, -1e9, 1e9, 0.7);
      CustomizeHistogramsYRange({hCutVarResidual}, false, -1e9, 1e9, 0.7);
      if(isMc) {
        ScalePlotVertically(hMcRatio, histoMc, 2);
        hMcRatio->GetYaxis()->SetTitle("MC / Fit");
        DivideHistoByFunction(hMcRatio, fitMc, integralOption);
        ScalePlotVertically(hMcResidual, histoMc, 2);
        hMcResidual->GetYaxis()->SetTitle("(MC - Fit) / #sigma_{MC}");
        EvalNormDifferenceHistoFromFunction(hMcResidual, fitMc, integralOption);
        CustomizeHistogramsYRange({hCutVarRatio, hMcRatio});
        CustomizeHistogramsYRange({hCutVarResidual, hMcResidual});
      }

      TCanvas ccRatio("ccRatio", "");
      ScaleCanvasVertically(&ccRatio, &ccFit, 2);
      TF1 oneline("oneline", "[0]", 0, 2);
      oneline.SetParameter(0, 1);
      oneline.SetLineColor(kBlack);
      oneline.SetLineStyle(7);
      if(isMc) hMcRatio->Draw("HIST PE");
      hCutVarRatio->Draw("HIST PE same");
      oneline.Draw("same");
      ccRatio.Print("ctfit.pdf", "pdf");

      TCanvas ccResidual("ccResidual", "");
      ScaleCanvasVertically(&ccResidual, &ccFit, 2);
      TF1 zeroline("oneline", "[0]", 0, 2);
      zeroline.SetParameter(0, 0);
      zeroline.SetLineColor(kBlack);
      zeroline.SetLineStyle(7);
      if(isMc) hMcResidual->Draw("HIST PE");
      hCutVarResidual->Draw("HIST PE same");
      zeroline.Draw("same");
      ccResidual.Print("ctfit.pdf", "pdf");

      if(IsSaveCanvasAsRoot && dropSet == 0) {
        fileOut->cd();
        ccRatio.Write();
        ccResidual.Write();
      }
    };

    for(int iDropSet=0; iDropSet<nDropSets; ++iDropSet) {
      const auto binsToDrop = HelperCtFit::EvalBinsToDrop(iDropSet);
      const auto binsToDropSize = static_cast<int>(binsToDrop.size());
      if(RunMode == HelperCtFit::RunModes::AllPossible && binsToDropSize >= nBins-1) continue;
      if(RunMode == HelperCtFit::RunModes::AllWoOne && binsToDropSize > 1) continue;
      TH1* hReco = dynamic_cast<TH1*>(hCutVarDiff->Clone());
      TH1* hMc = isMc ? dynamic_cast<TH1*>(hMcDiff->Clone()) : nullptr;

      for(const auto& binToDrop : binsToDrop) {
        HelperCtFit::ExcludeBin(hReco, binToDrop);
        if(isMc) HelperCtFit::ExcludeBin(hMc, binToDrop);
      }

      FitLifetimeHistos(hReco, hMc, iDropSet);
    } // nDropSets
    emptycanvas.Print("ctfit.pdf]", "pdf");

    if(IsSaveCanvasAsRoot) fileOut->Close();

    emptycanvas.Print("dropSummary.pdf[", "pdf");
    for(const auto& gr : graphsQa) {
      HelperCtFit::PrintCanvasDrops(gr, nBins);
    }
    emptycanvas.Print("dropSummary.pdf]", "pdf");

  } // promptnesses
}

int main(int argc, char* argv[]) {

  using namespace boost::program_options;
  options_description desc ("Allowed options");
  desc.add_options()
    ("help,h", "Print usage message")
    ("input-file,i", value<std::string>()->required(), "file with cut variation output")
    ("mc-file,m", value<std::string>()->default_value(""), "file with MC")
    ("syst-errors-file,s", value<std::string>()->default_value(""), "file with syst errors")
  ;

  variables_map args;
  store(parse_command_line(argc, argv, desc), args);

  if (args.count("help") || argc < 2) {
    std::cout << desc << "\n";
    return 0;
  }
  notify (args);

  const std::string fileNameCutVar = args["input-file"].as<std::string>();
  const std::string fileNameMC = args["mc-file"].as<std::string>();
  const std::string fileNameSystUnc = args["syst-errors-file"].as<std::string>();

  corrected_yields_qa2(fileNameCutVar, fileNameMC, fileNameSystUnc);

  return 0;
}
