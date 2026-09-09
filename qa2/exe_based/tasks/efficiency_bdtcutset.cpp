//
// Created by oleksii on 13.06.25.
//
#include "HelperGeneral.hpp"
#include "HelperMath.hpp"
#include "HelperPlot.hpp"

#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TString.h>
#include <TStyle.h>

#include <iostream>
#include <vector>

using namespace HelperGeneral;
using namespace HelperMath;
using namespace HelperPlot;

void RebinHistoToEdges(TH1*& histo, const std::vector<double>& edges);
void RebinHistoToEdges(TH2*& histo, const std::vector<double>& edges);

void efficiency_bdtcutset(const std::string& fileName) {
  LoadMacro("styles/mc_qa2.style.cc");
  gStyle->SetMarkerSize(1);
  gStyle->SetLineWidth(1);

  TFile* fileIn = OpenFileWithNullptrCheck(fileName);

  const std::string fileOutName = "efficiency_summary";

  // ========================= Configuration =================================
  const std::vector<double> lifeTimeRanges = {0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8};
  const std::vector<double> pTRanges = {1, 2, 3, 4, 5, 8, 12, 20};

  std::vector<float> bdtScores{0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90};
//   for(int i=0; i<=99; i++) {
//     bdtScores.emplace_back(0.01 * i);
//   }
// //   bdtScores.emplace_back(0.01);

  std::vector<std::pair<double, double>> pTIntervals{};
  for(size_t iPt = 0, nPts = pTRanges.size() - 1; iPt < nPts; ++iPt) {
//     pTIntervals.emplace_back(std::make_pair(pTRanges.at(iPt), pTRanges.at(iPt+1)));
  }

//   pTIntervals.emplace_back(std::make_pair(1., 20.));
  pTIntervals.emplace_back(std::make_pair(3., 20.));
//   pTIntervals.emplace_back(std::make_pair(4., 20.));

  // ==========================================================================

  const std::string tarSigShortcut = "NP";
  const std::vector<std::string> promptnesses{"prompt", "nonprompt"};
  const std::vector<std::string> weightPresences{"", "_W"};

  TFile* fileOut = TFile::Open((fileOutName + ".root").c_str(), "recreate");
  HelperMath::tensor<TGraphErrors*, 4> grEff = make_tensor<TGraphErrors*, 4>({promptnesses.size(), weightPresences.size(), pTIntervals.size(), lifeTimeRanges.size()-1}, nullptr);

  auto PtRangeString = [] (const std::pair<double, double>& pTInterval) {
    return "pT_" + HelperGeneral::to_string_with_precision(pTInterval.first, 0) + "_" + HelperGeneral::to_string_with_precision(pTInterval.second, 0);
  };

  auto PtRangeTitle = [] (const std::pair<double, double>& pTInterval) {
    return "#it{p}_{T}#in (" +
           HelperGeneral::to_string_with_precision(pTInterval.first, 0) + "#; " +
           HelperGeneral::to_string_with_precision(pTInterval.second, 0) + ") GeV/#it{c}";
  };

  for(size_t iPromptness=0, nPromptnesses=promptnesses.size(); iPromptness<nPromptnesses; ++iPromptness) {
    const std::string& promptness = promptnesses.at(iPromptness);
    std::cout << "Processing " << promptness << "\n";
    for(size_t iWeightPresence=0, nWeightPresences=weightPresences.size(); iWeightPresence < nWeightPresences; ++iWeightPresence) {
      for (size_t iPt = 0, nPts = pTIntervals.size(); iPt < nPts; ++iPt) {
        std::cout << "Processing iPt = " << iPt << "\t" << PtRangeString(pTIntervals.at(iPt)) << "\n";
        TH1* histoGen = GetObjectWithNullptrCheck<TH1>(fileIn, "gen/" + promptness + "/" + PtRangeString(pTIntervals.at(iPt)) + "/hT" + weightPresences.at(iWeightPresence));
        RebinHistoToEdges(histoGen, lifeTimeRanges);
        histoGen->UseCurrentStyle();
        CD(fileOut, "yields/" + promptness + "/" + PtRangeString(pTIntervals.at(iPt)));
        histoGen->Write(("gen" + weightPresences.at(iWeightPresence)).c_str());

        for (const auto& score: bdtScores) {
          const std::string sScore = to_string_with_precision(score, 2);
          if(score == bdtScores.at(0)) std::cout << "Processing iLifeTimeRange ";
          for (size_t iLifeTimeRange = 0; iLifeTimeRange < lifeTimeRanges.size() - 1 && score == bdtScores.at(0); ++iLifeTimeRange) {
            std::cout << iLifeTimeRange << " ";
            const Color_t grColor = promptness == "prompt" ? kRed : kBlue;
            const Style_t grLineStyle = weightPresences.at(iWeightPresence) == "" ? 1 : 7;
            const Style_t grMarkerStyle = weightPresences.at(iWeightPresence) == "" ? kFullSquare : kOpenSquare;
            grEff.at(iPromptness).at(iWeightPresence).at(iPt).at(iLifeTimeRange) = new TGraphErrors();
            auto gr = grEff.at(iPromptness).at(iWeightPresence).at(iPt).at(iLifeTimeRange);
            gr->SetName(("grEff_" + promptness + "_vs_" + tarSigShortcut + "_T" +
                         std::to_string(iLifeTimeRange)).c_str());
            gr->SetTitle(("bin #" + std::to_string(iLifeTimeRange + 1) + "#; T#in (" +
                          to_string_with_precision(lifeTimeRanges.at(iLifeTimeRange), 2) + "#; " +
                          to_string_with_precision(lifeTimeRanges.at(iLifeTimeRange + 1), 2) + ") ps#; " +
                          PtRangeTitle(pTIntervals.at(iPt))).c_str());
            gr->SetMarkerColor(grColor);
            gr->SetLineColor(grColor);
            gr->SetLineStyle(grLineStyle);
            gr->SetMarkerStyle(grMarkerStyle);
            gr->GetXaxis()->SetTitle(("bdt score " + tarSigShortcut).c_str());
            gr->GetYaxis()->SetTitle("Eff. #times Acc.");
            gr->SetMinimum(1e-5);
            gr->SetMaximum(1);
          }
          if(score == bdtScores.at(0)) std::cout << "\n";

          const std::string histoRecName{"rec/" + promptness + "/" + PtRangeString(pTIntervals.at(iPt)) + "/hT_" + tarSigShortcut + "gt" + sScore + weightPresences.at(iWeightPresence)};
          std::string histoRecName2D{histoRecName};
          ReplaceSubstrInStr(histoRecName2D, "hT", "h2T");
          TH1* histoRec = GetObjectWithNullptrCheck<TH1>(fileIn, histoRecName);
          TH2* histoRec2D = GetObjectWithNullptrCheck<TH2>(fileIn, histoRecName2D);
          RebinHistoToEdges(histoRec, lifeTimeRanges);
          RebinHistoToEdges(histoRec2D, lifeTimeRanges);
          histoRec->UseCurrentStyle();
          histoRec2D->UseCurrentStyle();

          CD(fileOut, "yields/" + promptness + "/" + PtRangeString(pTIntervals.at(iPt)));
          histoRec->Write(("rec_" + tarSigShortcut + "gt" + sScore + weightPresences.at(iWeightPresence)).c_str());
          histoRec2D->Write(("rec2D_" + tarSigShortcut + "gt" + sScore + weightPresences.at(iWeightPresence)).c_str()); // NOTE to start with

          const auto [histoEff, histoEffRelErr] = EvaluateEfficiencyHisto(histoRec, histoGen);
          const auto [histoResp, histoRespRelErr] = EvaluateResponseMatrix(histoRec2D, histoGen);

          CD(fileOut, "effs/" + promptness + "/" + PtRangeString(pTIntervals.at(iPt)));
          histoEff->Write(("eff_" + tarSigShortcut + "gt" + sScore + weightPresences.at(iWeightPresence)).c_str());
          histoResp->Write(("r_" + tarSigShortcut + "gt" + sScore + weightPresences.at(iWeightPresence)).c_str());

          CD(fileOut, "errs/" + promptness + "/" + PtRangeString(pTIntervals.at(iPt)));
          histoEffRelErr->Write(("err_" + tarSigShortcut + "gt" + sScore + weightPresences.at(iWeightPresence)).c_str());
          histoRespRelErr->Write(("r_err_" + tarSigShortcut + "gt" + sScore + weightPresences.at(iWeightPresence)).c_str());

          const std::string openOption = iPromptness + iPt + iWeightPresence == 0 ? "recreate" : "update";
          TFile* fileOutScore = TFile::Open(("Eff_times_Acc_Lc." + tarSigShortcut + "gt" + sScore + ".root").c_str(), openOption.c_str());
          CD(fileOutScore, PtRangeString(pTIntervals.at(iPt)));
          histoEff->Write((promptness + weightPresences.at(iWeightPresence)).c_str());
          for (size_t iLifeTimeRange = 0; iLifeTimeRange < lifeTimeRanges.size() - 1; ++iLifeTimeRange) {
            auto gr = grEff.at(iPromptness).at(iWeightPresence).at(iPt).at(iLifeTimeRange);
            gr->SetPoint(gr->GetN(), score, histoEff->GetBinContent(iLifeTimeRange + 1));
            gr->SetPointError(gr->GetN() - 1, 0, histoEff->GetBinError(iLifeTimeRange + 1));
          }
          fileOutScore->Close();
        } // bdtSignalLowerValues
      } // pTRanges
    } // weightPresences
  } // promptnesses

  TLegend leg(0.7, 0.7, 0.9, 0.9);
  size_t iPromptness = 0;
  for(const auto& promptness : promptnesses) {
    size_t iWeightPresence = 0;
    for(const auto& weightPresence : weightPresences) {
      leg.AddEntry(grEff.at(iPromptness).at(iWeightPresence).at(0).at(0), (promptness + weightPresence).c_str(), "PL");
      ++iWeightPresence;
    }
    ++iPromptness;
  }

  for(size_t iPt=0, nPts=pTIntervals.size(); iPt<nPts; ++iPt) {
    for (size_t iLifeTimeRange = 0; iLifeTimeRange < lifeTimeRanges.size() - 1; ++iLifeTimeRange) {
      const std::string priBra = lifeTimeRanges.size() - 1 == 1 ? "" : iLifeTimeRange == 0 ? "(" : iLifeTimeRange == lifeTimeRanges.size() - 2 ? ")" : "";
      TCanvas cc("cc", "");
      cc.SetCanvasSize(1200, 800);
      cc.SetLogy();
      grEff.at(0).at(0).at(iPt).at(iLifeTimeRange)->Draw("APE"); // 0 stands for prompt, 0 for no-weights
      grEff.at(0).at(1).at(iPt).at(iLifeTimeRange)->Draw("PE same"); // 0 stands for prompt, 1 for weights
      grEff.at(1).at(0).at(iPt).at(iLifeTimeRange)->Draw("PE same"); // 1 stands for nonprompt, 0 for no-weights
      grEff.at(1).at(1).at(iPt).at(iLifeTimeRange)->Draw("PE same"); // 1 stands for nonprompt, 1 for weights
      leg.Draw("same");
      cc.Print(("grEff_vs_" + tarSigShortcut + "_" + PtRangeString(pTIntervals.at(iPt)) + ".pdf" + priBra).c_str());
    } // lifeTimeRanges
  } // pTRanges
  fileOut->Close();
  fileIn->Close();
}

void RebinHistoToEdges(TH1*& histo, const std::vector<double>& edges) { // TODO consider mv to Helper
  histo = dynamic_cast<TH1*>(histo->Rebin(edges.size() - 1, histo->GetName(), edges.data()));
}

void RebinHistoToEdges(TH2*& histo, const std::vector<double>& edges) { // TODO consider mv to Helper
  const int nBins = edges.size() - 1;

  // Preserve whether Sumw2 was actually enabled.
  const bool hasSumw2 = histo->GetSumw2N() > 0;

  auto* rebinned = new TH2D("", histo->GetTitle(), nBins, edges.data(), nBins, edges.data());
  rebinned->SetDirectory(nullptr);
  rebinned->SetName(histo->GetName());

  rebinned->GetXaxis()->SetTitle(histo->GetXaxis()->GetTitle());
  rebinned->GetYaxis()->SetTitle(histo->GetYaxis()->GetTitle());

  if (hasSumw2) rebinned->Sumw2();

  // Include underflow and overflow, like histogram rebinning should.
  for (int ix = 0; ix <= histo->GetNbinsX() + 1; ++ix) {
    const double x = histo->GetXaxis()->GetBinCenter(ix);

    int jx;
    if (ix == 0) jx = 0;
    else if (ix == histo->GetNbinsX() + 1) jx = nBins + 1;
    else jx = rebinned->GetXaxis()->FindBin(x);

    for (int iy = 0; iy <= histo->GetNbinsY() + 1; ++iy) {
      const double y = histo->GetYaxis()->GetBinCenter(iy);

      int jy;
      if (iy == 0) jy = 0;
      else if (iy == histo->GetNbinsY() + 1) jy = nBins + 1;
      else jy = rebinned->GetYaxis()->FindBin(y);

      const int oldBin = histo->GetBin(ix, iy);
      const int newBin = rebinned->GetBin(jx, jy);

      rebinned->SetBinContent(newBin, rebinned->GetBinContent(newBin) + histo->GetBinContent(oldBin));

      if (hasSumw2) {
        const double oldErr = histo->GetBinError(oldBin);
        const double newErr = rebinned->GetBinError(newBin);

        rebinned->SetBinError(newBin, std::sqrt(newErr * newErr + oldErr * oldErr));
      }
    } // iy
  } // ix

  // Preserve number of entries.
  rebinned->SetEntries(histo->GetEntries());

  delete histo;
  histo = rebinned;
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
