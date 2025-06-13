//
// Created by oleksii on 13.06.25.
//
#include "Helper.hpp"
#include "HelperEfficiency.hpp"

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TString.h>

#include <iostream>

using namespace Helper;
using namespace HelperEfficiency;

void RebinHistoToEdges(TH1*& histo, const std::vector<double>& edges);

void efficiency_bdtdiff(const std::string& fileName, bool isSaveToRoot) {
  LoadMacro("styles/mc_qa2.style.cc");

  TFile* fileIn = OpenFileWithNullptrCheck(fileName);

  const std::string fileOutName = "efficiency_bdtdiff";

  // ========================= Configuration =================================
  const std::string signalShortcut = "P";
  const std::vector<double> lifeTimeRanges = {0.2, 0.35, 0.5, 0.7, 0.9, 1.6};

  std::vector<float> bdtSignalLowerValues;
  for(int iB=0; iB<=50; iB++) {
    bdtSignalLowerValues.emplace_back(0.02 * iB);
  }
  // ==========================================================================

  const std::vector<std::string> promptnesses{"prompt", "nonprompt"};

  TFile* fileOut{nullptr};
  if(isSaveToRoot) fileOut = TFile::Open((fileOutName + ".root").c_str(), "recreate");

  for(const auto& promptness : promptnesses) {
    TH1* histoGen = GetObjectWithNullptrCheck<TH1>(fileIn, "gen/" + promptness + "/hT");
    RebinHistoToEdges(histoGen, lifeTimeRanges);
    histoGen->UseCurrentStyle();
    TCanvas cGen("cGen", "");
    cGen.SetCanvasSize(1200, 800);
    histoGen->Draw();
    AddOneLineText(promptness, {0.35, 0.95, 0.45, 0.99});
    cGen.Print(("yieldGen." + promptness + ".pdf").c_str(), "pdf");

    if(isSaveToRoot) {
      CD(fileOut, "yields/" + promptness);
      histoGen->Write("gen");
    }

    std::string priBra = "(";
    for(const auto& score : bdtSignalLowerValues) {
      const std::string sScore = to_string_with_precision(score, 2);

      TH1* histoRec = GetObjectWithNullptrCheck<TH1>(fileIn, "rec/" + promptness + "/hT_" + signalShortcut + "gt" + sScore);
      RebinHistoToEdges(histoRec, lifeTimeRanges);
      histoRec->UseCurrentStyle();
      TCanvas cRec("cRec", "");
      cRec.SetCanvasSize(1200, 800);
      histoRec->Draw();
      AddOneLineText(promptness, {0.35, 0.95, 0.45, 0.99});
      AddOneLineText("bdt(" + signalShortcut + ") > " + sScore, {0.65, 0.95, 0.75, 0.99});
      cRec.Print(("yieldRec." + promptness + ".pdf" + priBra).c_str(), "pdf");

      auto [histoEff, histoEffRelErr] = EvaluateEfficiencyHisto(histoRec, histoGen);

      TCanvas cEff("cEff", "");
      cEff.SetCanvasSize(1200, 800);
      histoEff->Draw();
      AddOneLineText(promptness, {0.35, 0.95, 0.45, 0.99});
      AddOneLineText("bdt(" + signalShortcut + ") > " + sScore, {0.65, 0.95, 0.75, 0.99});
      cEff.Print(("eff." + promptness + ".pdf" + priBra).c_str(), "pdf");

      TCanvas cErr("cErr", "");
      cErr.SetCanvasSize(1200, 800);
      histoEffRelErr->Draw();
      AddOneLineText(promptness, {0.35, 0.95, 0.45, 0.99});
      AddOneLineText("bdt(" + signalShortcut + ") > " + sScore, {0.65, 0.95, 0.75, 0.99});
      cErr.Print(("err." + promptness + ".pdf" + priBra).c_str(), "pdf");

      if(isSaveToRoot) {
        CD(fileOut, "yields/" + promptness);
        histoRec->Write(("rec_" + signalShortcut + "gt" + sScore).c_str());
        CD(fileOut, "effs/" + promptness);
        histoEff->Write(("eff_" + signalShortcut + "gt" + sScore).c_str());
        CD(fileOut, "errs/" + promptness);
        histoEffRelErr->Write(("err_" + signalShortcut + "gt" + sScore).c_str());
      }

      priBra = "";
    } // bdtSignalLowerValues

    CloseCanvasPrinting({"yieldRec." + promptness, "eff." + promptness, "err." + promptness});

  } // promptnesses

  if(isSaveToRoot) fileOut->Close();
  fileIn->Close();
}

void RebinHistoToEdges(TH1*& histo, const std::vector<double>& edges) {
  histo = dynamic_cast<TH1*>(histo->Rebin(edges.size() - 1,histo->GetName(),edges.data()));
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./efficiency_bdtdiff fileName (isSaveRoot=false)" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileName = argv[1];
  const bool isSaveToRoot = argc > 2 ? string_to_bool(argv[2]) : false;

  efficiency_bdtdiff(fileName, isSaveToRoot);

  return 0;
}
