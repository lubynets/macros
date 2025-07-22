//
// Created by oleksii on 13.06.25.
//
#include "HelperGeneral.hpp"
#include "HelperMath.hpp"
#include "HelperPlot.hpp"

#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TLegend.h>
#include <TString.h>

#include <iostream>

using namespace HelperGeneral;
using namespace HelperMath;
using namespace HelperPlot;

void RebinHistoToEdges(TH1*& histo, const std::vector<double>& edges);

void efficiency_bdtcutset(const std::string& fileName) {
  LoadMacro("styles/mc_qa2.style.cc");

  TFile* fileIn = OpenFileWithNullptrCheck(fileName);

  const std::string fileOutName = "efficiency_summary";

  // ========================= Configuration =================================
  const std::vector<double> lifeTimeRanges = {0.2, 0.35, 0.5, 0.7, 0.9, 1.6};
  const std::vector<double> pTRanges = {0, 2, 5, 8, 12, 20};

  std::vector<float> bdtScores;
  for(int i=0; i<=99; i++) {
    bdtScores.emplace_back(0.01 * i);
  }

  const std::string tarSigShortcut = "NP";
  // ==========================================================================

  const std::vector<std::string> promptnesses{"prompt", "nonprompt"};

  TFile* fileOut = TFile::Open((fileOutName + ".root").c_str(), "recreate");
  HelperMath::tensor<TGraphErrors*, 3> grEff = make_tensor<TGraphErrors*, 3>({promptnesses.size(), pTRanges.size(), lifeTimeRanges.size()-1}, nullptr);

  auto PtRangeString = [&] (size_t iPt) {
    std::pair<size_t, size_t> iPTMinMax = (iPt == pTRanges.size()-1) ? std::pair<size_t, size_t>{0, pTRanges.size()-1} : std::pair<size_t, size_t>{iPt, iPt+1};
    return "pT_" + HelperGeneral::to_string_with_precision(pTRanges.at(iPTMinMax.first), 0) + "_" + HelperGeneral::to_string_with_precision(pTRanges.at(iPTMinMax.second), 0);
  };

  auto PtRangeTitle = [&] (size_t iPt) {
    std::pair<size_t, size_t> iPTMinMax = (iPt == pTRanges.size()-1) ? std::pair<size_t, size_t>{0, pTRanges.size()-1} : std::pair<size_t, size_t>{iPt, iPt+1};
    return "#it{p}_{T}#in (" +
           HelperGeneral::to_string_with_precision(pTRanges.at(iPTMinMax.first), 0) + "#; " +
           HelperGeneral::to_string_with_precision(pTRanges.at(iPTMinMax.second), 0) + ") GeV/#it{c}";
  };

  for(size_t iPromptness=0, nPromptnesses=promptnesses.size(); iPromptness<nPromptnesses; ++iPromptness) {
    const std::string& promptness = promptnesses.at(iPromptness);
    std::cout << "Processing " << promptness << "\n";
    for(size_t iPt=0, nPts=pTRanges.size()-1; iPt<nPts; ++iPt) {
      std::cout << "Processing iPt = " << iPt << "\t" << PtRangeString(iPt) << "\n";
      TH1* histoGen = GetObjectWithNullptrCheck<TH1>(fileIn, "gen/" + promptness + "/" + PtRangeString(iPt) + "/hT");
      RebinHistoToEdges(histoGen, lifeTimeRanges);
      histoGen->UseCurrentStyle();
      CD(fileOut, "yields/" + promptness + "/" + PtRangeString(iPt));
      histoGen->Write("gen");

        for (const auto& score : bdtScores) {
          const std::string sScore = to_string_with_precision(score, 2);
          for (int iLifeTimeRange = 0; iLifeTimeRange < lifeTimeRanges.size() - 1 && score == bdtScores.at(0); ++iLifeTimeRange) {
            std::cout << "Processing iLifeTimeRange " << iLifeTimeRange << "\n";
            const Color_t grColor = promptness == "prompt" ? kRed : kBlue;
            grEff.at(iPromptness).at(iPt).at(iLifeTimeRange) = new TGraphErrors();
            auto gr = grEff.at(iPromptness).at(iPt).at(iLifeTimeRange);
            gr->SetName(("grEff_" + promptness + "_vs_" + tarSigShortcut + "_T" +
                         std::to_string(iLifeTimeRange)).c_str());
            gr->SetTitle(("bin #" + std::to_string(iLifeTimeRange + 1) + "#; T#in (" +
                          to_string_with_precision(lifeTimeRanges.at(iLifeTimeRange), 2) + "#; " +
                          to_string_with_precision(lifeTimeRanges.at(iLifeTimeRange + 1), 2) + ") ps#; " +
                          PtRangeTitle(iPt)).c_str());
            std::cout << PtRangeTitle(iPt) << "\t";
            std::cout << gr->GetTitle() << "\n";
            gr->SetMarkerColor(grColor);
            gr->SetLineColor(grColor);
            gr->GetXaxis()->SetTitle(("bdt score " + tarSigShortcut).c_str());
            gr->GetYaxis()->SetTitle("Eff. #times Acc.");
            gr->SetMinimum(1e-5);
            gr->SetMaximum(1);
          }

          TH1* histoRec = GetObjectWithNullptrCheck<TH1>(fileIn,"rec/" + promptness + "/" + PtRangeString(iPt) + "/hT_" + tarSigShortcut + "gt" + sScore);
          RebinHistoToEdges(histoRec, lifeTimeRanges);
          histoRec->UseCurrentStyle();

          CD(fileOut, "yields/" + promptness + "/" + PtRangeString(iPt));
          histoRec->Write(("rec_" + tarSigShortcut + "gt" + sScore).c_str());

          auto [histoEff, histoEffRelErr] = EvaluateEfficiencyHisto(histoRec, histoGen);

          CD(fileOut, "effs/" + promptness + "/" + PtRangeString(iPt));
          histoEff->Write(("eff_" + tarSigShortcut + "gt" + sScore).c_str());

          CD(fileOut, "errs/" + promptness + "/" + PtRangeString(iPt));
          histoEffRelErr->Write(("err_" + tarSigShortcut + "gt" + sScore).c_str());

          const std::string openOption = promptness == "prompt" ? "recreate" : "update";
          TFile* fileOutScore = TFile::Open(("Eff_times_Acc_Lc." + tarSigShortcut + "gt" + sScore + ".root").c_str(),openOption.c_str());
          CD(fileOutScore, PtRangeString(iPt));
          histoEff->Write(promptness.c_str());
          for (int iLifeTimeRange = 0; iLifeTimeRange < lifeTimeRanges.size() - 1; ++iLifeTimeRange) {
            auto gr = grEff.at(iPromptness).at(iPt).at(iLifeTimeRange);
            gr->SetPoint(gr->GetN(), score, histoEff->GetBinContent(iLifeTimeRange + 1));
            gr->SetPointError(gr->GetN() - 1, 0, histoEff->GetBinError(iLifeTimeRange + 1));
          }
          fileOutScore->Close();
        } // bdtSignalLowerValues
    } // pTRanges
  } // promptnesses

  TLegend leg(0.7, 0.7, 0.9, 0.9);
  size_t iPromptness = 0;
  for(const auto& promptness : promptnesses) {
    leg.AddEntry(grEff.at(iPromptness).at(0).at(0), promptness.c_str(), "PL");
    ++iPromptness;
  }

  for(size_t iPt=0, nPts=pTRanges.size()-1; iPt<nPts; ++iPt) {
    for (int iLifeTimeRange = 0; iLifeTimeRange < lifeTimeRanges.size() - 1; ++iLifeTimeRange) {
      const std::string priBra = lifeTimeRanges.size() - 1 == 1 ? "" : iLifeTimeRange == 0 ? "(" : iLifeTimeRange == lifeTimeRanges.size() - 2 ? ")" : "";
      TCanvas cc("cc", "");
      cc.SetCanvasSize(1200, 800);
      cc.SetLogy();
      grEff.at(0).at(iPt).at(iLifeTimeRange)->Draw("APE"); // 0 stands for prompt
      grEff.at(1).at(iPt).at(iLifeTimeRange)->Draw("PE same"); // 1 stands for nonprompt
      leg.Draw("same");
      cc.Print(("grEff_vs_" + tarSigShortcut + "_" + PtRangeString(iPt) + ".pdf" + priBra).c_str());
    } // lifeTimeRanges
  } // pTRanges
  fileOut->Close();
  fileIn->Close();
}

void RebinHistoToEdges(TH1*& histo, const std::vector<double>& edges) {
  histo = dynamic_cast<TH1*>(histo->Rebin(edges.size() - 1,histo->GetName(),edges.data()));
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./efficiency_bdtcutset fileName" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileName = argv[1];

  efficiency_bdtcutset(fileName);

  return 0;
}
