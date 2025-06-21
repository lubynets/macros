//
// Created by oleksii on 13.06.25.
//
#include "Helper.hpp"
#include "HelperEfficiency.hpp"

#include <TFile.h>
#include <TH1.h>
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
  for(int i=0; i<=19; i++) {
    targetSignalsFill.emplace_back(0.05*i);
  }

  const std::vector<TargetSignal> targetSignals {
    {"P", targetSignalsFill},
    {"NP", targetSignalsFill},
  };
  // ==========================================================================

  const std::vector<std::string> promptnesses{"prompt", "nonprompt"};

  TFile* fileOut = TFile::Open((fileOutName + ".root").c_str(), "recreate");

  bool isScoresPrinted{false};
  for(const auto& promptness : promptnesses) {
    TH1* histoGen = GetObjectWithNullptrCheck<TH1>(fileIn, "gen/" + promptness + "/hT");
    RebinHistoToEdges(histoGen, lifeTimeRanges);
    histoGen->UseCurrentStyle();
    CD(fileOut, "yields/" + promptness);
    histoGen->Write("gen");

    for(const auto& tarSig : targetSignals) {
      if(!isScoresPrinted) std::cout << "SCORES_" << tarSig.shortcut_ << "=\"";
      for (const auto& score: tarSig.bdt_lower_values_) {
        if(!isScoresPrinted) std::cout << to_string_with_precision(score, 2) << " ";
        const std::string sScore = to_string_with_precision(score, 2);

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
        fileOutScore->Close();
      } // bdtSignalLowerValues
      if(!isScoresPrinted) std::cout << "\"\n";
    } // targetSignals
    isScoresPrinted = true;
  } // promptnesses

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
