//
// Created by oleksii on 13.06.25.
//
#include "Helper.hpp"
#include "HelperEfficiency.hpp"

#include <TFile.h>
#include <TH1.h>
#include <TLegend.h>
#include <TString.h>

#include <iostream>

using namespace Helper;
using namespace HelperEfficiency;

void RebinHistoToEdges(TH1*& histo, const std::vector<double>& edges);
struct TargetSignal {
  std::string shortcut_;
  std::vector<float> bdt_lower_values_;
};

void efficiency_bdtcutset(const std::string& fileName) {
  LoadMacro("styles/mc_qa2.style.cc");

  TFile* fileIn = OpenFileWithNullptrCheck(fileName);

  const std::string fileOutName = "efficiency_summary";

  // ========================= Configuration =================================
  const std::vector<double> lifeTimeRanges = {0.2, 0.35, 0.5, 0.7, 0.9, 1.6};

  std::vector<float> targetSignalsFill;
  for(int i=0; i<=99; i++) {
    targetSignalsFill.emplace_back(0.01*i);
  }

  const std::vector<TargetSignal> targetSignals {
    {"P", targetSignalsFill},
    {"NP", targetSignalsFill},
  };
  // ==========================================================================

  const std::vector<std::string> promptnesses{"prompt", "nonprompt"};

  TFile* fileOut = TFile::Open((fileOutName + ".root").c_str(), "recreate");
  std::vector<std::vector<std::vector<TGraphErrors*>>> grEff(lifeTimeRanges.size()-1);
  for(auto& gr : grEff) {
    gr.resize(promptnesses.size());
    for(auto& grr : gr) {
      grr.resize(targetSignals.size());
    }
  }

  int iPromptness{0};
  for(const auto& promptness : promptnesses) {
    TH1* histoGen = GetObjectWithNullptrCheck<TH1>(fileIn, "gen/" + promptness + "/hT");
    RebinHistoToEdges(histoGen, lifeTimeRanges);
    histoGen->UseCurrentStyle();
    CD(fileOut, "yields/" + promptness);
    histoGen->Write("gen");

    int iTargetSignal{0};
    for(const auto& tarSig : targetSignals) {
      for (const auto& score: tarSig.bdt_lower_values_) {
        const std::string sScore = to_string_with_precision(score, 2);
        for(int iLifeTimeRange=0; iLifeTimeRange<lifeTimeRanges.size()-1 && score == tarSig.bdt_lower_values_.at(0); ++iLifeTimeRange){
          const Color_t grColor = promptness == "prompt" ? kRed : kBlue;
          grEff.at(iLifeTimeRange).at(iPromptness).at(iTargetSignal) = new TGraphErrors();
          auto gr = grEff.at(iLifeTimeRange).at(iPromptness).at(iTargetSignal);
          gr->SetName(("grEff_" + promptness + "_vs_" + tarSig.shortcut_ + "_T" + std::to_string(iLifeTimeRange)).c_str());
          gr->SetTitle(("bin #" + std::to_string(iLifeTimeRange) + "#; T#in (" + to_string_with_precision(lifeTimeRanges.at(iLifeTimeRange), 2) + "#; " + to_string_with_precision(lifeTimeRanges.at(iLifeTimeRange+1), 2) + ") ps").c_str());
          gr->SetMarkerColor(grColor);
          gr->SetLineColor(grColor);
          gr->GetXaxis()->SetTitle(("bdt score " + tarSig.shortcut_).c_str());
          gr->GetYaxis()->SetTitle("Eff. #times Acc.");
          gr->SetMinimum(1e-5);
          gr->SetMaximum(1);
        }

        TH1* histoRec = GetObjectWithNullptrCheck<TH1>(fileIn,"rec/" + promptness + "/hT_" + tarSig.shortcut_ + "gt" + sScore);
        RebinHistoToEdges(histoRec, lifeTimeRanges);
        histoRec->UseCurrentStyle();

        CD(fileOut, "yields/" + promptness);
        histoRec->Write(("rec_" + tarSig.shortcut_ + "gt" + sScore).c_str());

        auto [histoEff, histoEffRelErr] = EvaluateEfficiencyHisto(histoRec, histoGen);

        CD(fileOut, "effs/" + promptness);
        histoEff->Write(("eff_" + tarSig.shortcut_ + "gt" + sScore).c_str());

        CD(fileOut, "errs/" + promptness);
        histoEffRelErr->Write(("err_" + tarSig.shortcut_ + "gt" + sScore).c_str());

        const std::string openOption = promptness == "prompt" ? "recreate" : "update";
        TFile* fileOutScore = TFile::Open(("Eff_times_Acc_Lc." + tarSig.shortcut_ + "gt" + sScore + ".root").c_str(),
                                          openOption.c_str());
        histoEff->Write(promptness.c_str());
        for(int iLifeTimeRange=0; iLifeTimeRange<lifeTimeRanges.size()-1; ++iLifeTimeRange) {
          auto gr = grEff.at(iLifeTimeRange).at(iPromptness).at(iTargetSignal);
          gr->SetPoint(gr->GetN(), score, histoEff->GetBinContent(iLifeTimeRange+1));
          gr->SetPointError(gr->GetN()-1, 0, histoEff->GetBinError(iLifeTimeRange+1));
        }
        fileOutScore->Close();
      } // bdtSignalLowerValues
      ++iTargetSignal;
    } // targetSignals
    ++iPromptness;
  } // promptnesses

  TLegend leg(0.7, 0.7, 0.9, 0.9);
  iPromptness = 0;
  for(const auto& promptness : promptnesses) {
    leg.AddEntry(grEff.at(0).at(iPromptness).at(0), promptness.c_str(), "PL");
    ++iPromptness;
  }

  for(int iLifeTimeRange=0; iLifeTimeRange<lifeTimeRanges.size()-1; ++iLifeTimeRange) {
    const std::string priBra = lifeTimeRanges.size()-1 == 1 ? "" : iLifeTimeRange == 0 ? "(" : iLifeTimeRange == lifeTimeRanges.size()-2 ? ")" : "";
    for(int iTargetSignal=0; iTargetSignal<targetSignals.size(); ++iTargetSignal) {
      TCanvas cc("cc", "");
      cc.SetCanvasSize(1200, 800);
      cc.SetLogy();
      grEff.at(iLifeTimeRange).at(0).at(iTargetSignal)->Draw("APE");
      grEff.at(iLifeTimeRange).at(1).at(iTargetSignal)->Draw("PE same");
      leg.Draw("same");
      cc.Print(("grEff_vs_" + targetSignals.at(iTargetSignal).shortcut_ + ".pdf" + priBra).c_str());
    } // targetSignals
  } // lifeTimeRanges
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
