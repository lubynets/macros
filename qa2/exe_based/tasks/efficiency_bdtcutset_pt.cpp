//
// Created by oleksii on 25.11.26.
//
#include "HelperGeneral.hpp"
#include "HelperMath.hpp"
#include "HelperPlot.hpp"

#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>

#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <cstdlib>
#include <cmath>

using namespace HelperGeneral;
using namespace HelperMath;
using namespace HelperPlot;

std::string gCmdArgs{};
bool gIsSaveMeta{true};

void efficiency_bdtcutset(const std::string& fileName) {
  LoadMacro("styles/mc_qa2.style.cc");
  gStyle->SetMarkerSize(1);
  gStyle->SetLineWidth(1);

  TFile* fileIn = OpenFileWithNullptrCheck(fileName);
  const std::string fileOutName = "efficiency_summary";
  const std::string tarSigShortcut = "NP";
  const std::vector<std::string> promptnesses{"prompt", "nonprompt"};
  std::vector<std::string> weightPresences{""};

  const bool isPromptWeighted = fileIn->Get<TH1>("gen/prompt/hPt_W") != nullptr;
  const bool isNonPromptWeighted = fileIn->Get<TH1>("gen/nonprompt/hPt_W") != nullptr;
  if(isPromptWeighted != isNonPromptWeighted) throw std::runtime_error("Inconsistent weighted generated histograms");
  if(isPromptWeighted) weightPresences.emplace_back("_W");

  // Take the pT binning and the available BDT scores from the input file
  TH1* histoPt = GetObjectWithNullptrCheck<TH1>(fileIn, "gen/prompt/hPt");
  const size_t nPts = histoPt->GetNbinsX();
  if(nPts == 0) throw std::runtime_error("No pT bins found");
  std::vector<double> pTRanges{};
  for(size_t iPt=0; iPt<=nPts; ++iPt) {
    pTRanges.emplace_back(histoPt->GetXaxis()->GetBinLowEdge(iPt+1));
  }

  std::vector<double> bdtScores{};
  for(int i=0; i<=99; i++) {
    const double score = 0.01 * i;
    const std::string sScore = to_string_with_precision(score, 2);
    const std::string histoName = "rec/prompt/hPt_NPgt" + sScore;
    if(fileIn->Get<TH1>(histoName.c_str())) bdtScores.emplace_back(score);
  }
  if(bdtScores.empty()) throw std::runtime_error("No reconstructed BDT scan histograms found");

  auto CheckPtBinning = [&](const TH1* histo) {
    if(histo->GetDimension() != 1 || histo->GetNbinsX() != static_cast<int>(nPts)) throw std::runtime_error("Inconsistent pT binning");
    for(size_t iPt=0; iPt<=nPts; ++iPt) {
      if(std::abs(histo->GetXaxis()->GetBinLowEdge(iPt+1) - pTRanges.at(iPt)) > 1e-8) throw std::runtime_error("Inconsistent pT bin edges");
    }
  };

  auto PtRangeString = [&](size_t iPt) {
    return "pT_" + to_string_with_precision(pTRanges.at(iPt), 0) + "_" + to_string_with_precision(pTRanges.at(iPt+1), 0);
  };

  auto PtRangeTitle = [&](size_t iPt) {
    return "#it{p}_{T}#in (" +
           to_string_with_precision(pTRanges.at(iPt), 0) + "#; " +
           to_string_with_precision(pTRanges.at(iPt+1), 0) + ") GeV/#it{c}";
  };

  TFile* fileOut = TFile::Open((fileOutName + ".root").c_str(), "recreate");
  HelperMath::tensor<TGraphErrors*, 3> grEff = make_tensor<TGraphErrors*, 3>({promptnesses.size(), weightPresences.size(), nPts}, nullptr);

  for(size_t iPromptness=0, nPromptnesses=promptnesses.size(); iPromptness<nPromptnesses; ++iPromptness) {
    const std::string& promptness = promptnesses.at(iPromptness);
    std::cout << "Processing " << promptness << "\n";
    for(size_t iWeightPresence=0, nWeightPresences=weightPresences.size(); iWeightPresence < nWeightPresences; ++iWeightPresence) {
      const std::string& weightPresence = weightPresences.at(iWeightPresence);
      TH1* histoGen = GetObjectWithNullptrCheck<TH1>(fileIn, "gen/" + promptness + "/hPt" + weightPresence);
      CheckPtBinning(histoGen);
      histoGen->UseCurrentStyle();
      CD(fileOut, "yields/" + promptness);
      histoGen->Write(("gen" + weightPresence).c_str());

      for(size_t iPt=0; iPt<nPts; ++iPt) {
        const int grColor = promptness == "prompt" ? kRed : kBlue;
        const int grLineStyle = weightPresence == "" ? 1 : 7;
        const int grMarkerStyle = weightPresence == "" ? kFullSquare : kOpenSquare;
        grEff.at(iPromptness).at(iWeightPresence).at(iPt) = new TGraphErrors();
        auto gr = grEff.at(iPromptness).at(iWeightPresence).at(iPt);
        gr->SetName(("grEff_" + promptness + "_vs_" + tarSigShortcut + "_" + PtRangeString(iPt) + weightPresence).c_str());
        gr->SetTitle(PtRangeTitle(iPt).c_str());
        gr->SetMarkerColor(grColor);
        gr->SetLineColor(grColor);
        gr->SetLineStyle(grLineStyle);
        gr->SetMarkerStyle(grMarkerStyle);
        gr->GetXaxis()->SetTitle(("bdt score " + tarSigShortcut).c_str());
        gr->GetYaxis()->SetTitle("Eff. #times Acc.");
        gr->SetMinimum(1e-5);
        gr->SetMaximum(1);
      }

      for(const auto& score : bdtScores) {
        const std::string sScore = to_string_with_precision(score, 2);
        std::cout << "Processing bdt score = " << sScore << weightPresence << "\n";
        const std::string histoRecName{"rec/" + promptness + "/hPt_" + tarSigShortcut + "gt" + sScore + weightPresence};
        TH1* histoRec = GetObjectWithNullptrCheck<TH1>(fileIn, histoRecName);
        CheckPtBinning(histoRec);
        histoRec->UseCurrentStyle();

        CD(fileOut, "yields/" + promptness);
        histoRec->Write(("rec_" + tarSigShortcut + "gt" + sScore + weightPresence).c_str());

        const auto [histoEff, histoEffRelErr] = EvaluateEfficiencyHisto(histoRec, histoGen);
        histoEff->SetDirectory(nullptr);
        histoEffRelErr->SetDirectory(nullptr);
        CD(fileOut, "effs/" + promptness);
        histoEff->Write(("eff_" + tarSigShortcut + "gt" + sScore + weightPresence).c_str());
        CD(fileOut, "errs/" + promptness);
        histoEffRelErr->Write(("err_" + tarSigShortcut + "gt" + sScore + weightPresence).c_str());

        const std::string openOption = iPromptness + iWeightPresence == 0 ? "recreate" : "update";
        TFile* fileOutScore = TFile::Open(("Eff_times_Acc_Lc." + tarSigShortcut + "gt" + sScore + ".root").c_str(), openOption.c_str());
        fileOutScore->cd();
        histoEff->Write((promptness + weightPresence).c_str());
        for(size_t iPt=0; iPt<nPts; ++iPt) {
          auto gr = grEff.at(iPromptness).at(iWeightPresence).at(iPt);
          gr->SetPoint(gr->GetN(), score, histoEff->GetBinContent(iPt+1));
          gr->SetPointError(gr->GetN() - 1, 0, histoEff->GetBinError(iPt+1));
        }
        if (gIsSaveMeta && openOption == "recreate") HelperGeneral::EmbedMetaData(fileOutScore, "efficiency_bdtcutset_pt.cpp", gCmdArgs);
        fileOutScore->Close();
        delete fileOutScore;
        delete histoEff;
        delete histoEffRelErr;
      } // bdtScores
    } // weightPresences
  } // promptnesses

  TLegend leg(0.7, 0.7, 0.9, 0.9);
  for(size_t iPromptness=0; iPromptness<promptnesses.size(); ++iPromptness) {
    for(size_t iWeightPresence=0; iWeightPresence<weightPresences.size(); ++iWeightPresence) {
      leg.AddEntry(grEff.at(iPromptness).at(iWeightPresence).at(0), (promptnesses.at(iPromptness) + weightPresences.at(iWeightPresence)).c_str(), "PL");
    }
  }

  for(size_t iPt=0; iPt<nPts; ++iPt) {
    TCanvas cc("cc", "");
    cc.SetCanvasSize(1200, 800);
    cc.SetLogy();
    for(size_t iPromptness=0; iPromptness<promptnesses.size(); ++iPromptness) {
      for(size_t iWeightPresence=0; iWeightPresence<weightPresences.size(); ++iWeightPresence) {
        auto gr = grEff.at(iPromptness).at(iWeightPresence).at(iPt);
        gr->Draw(iPromptness + iWeightPresence == 0 ? "APE" : "PE same");
        CD(fileOut, "graphs/" + PtRangeString(iPt));
        gr->Write();
      }
    }
    leg.Draw("same");
    cc.Print(("grEff_vs_" + tarSigShortcut + "_" + PtRangeString(iPt) + ".pdf").c_str());
  } // pTRanges

  if (gIsSaveMeta) HelperGeneral::EmbedMetaData(fileOut, "efficiency_bdtcutset_pt.cpp", gCmdArgs);
  fileOut->Close();
  fileIn->Close();
  for(size_t iPromptness=0; iPromptness<promptnesses.size(); ++iPromptness) {
    for(size_t iWeightPresence=0; iWeightPresence<weightPresences.size(); ++iWeightPresence) {
      for(size_t iPt=0; iPt<nPts; ++iPt) {
        delete grEff.at(iPromptness).at(iWeightPresence).at(iPt);
      }
    }
  }
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./efficiency_bdtcutset_pt fileName" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileName = argv[1];

  for(int iArg=1; iArg<argc; ++iArg) {
    gCmdArgs += argv[iArg];
    gCmdArgs += " ";
  }

  efficiency_bdtcutset(fileName);

  return 0;
}
