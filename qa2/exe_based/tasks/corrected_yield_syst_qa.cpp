//
// Created by oleksii on 08.04.26.
//
#include "HelperGeneral.hpp"

#include <TFile.h>
#include <TH1.h>

#include <array>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

using namespace HelperGeneral;
using namespace std::string_literals;

template <typename T>
struct SystematicDof {
  std::string prefix_{};
  std::vector<T> values_{};
};

void corrected_yield_syst_qa() {
  LoadMacro("styles/mc_qa2.style.cc");

  const std::string meanCutVarFileName{"/lustre/alice/users/lubynets/syst/cutVar/outputs/CutVarLc.merged.root"};

  const std::string inputCommonPath{"/lustre/alice/users/lubynets/syst/cutVar/outputs"};
  SystematicDof<double> leftRanges{"lera", {2.12, 2.14}};
  SystematicDof<double> rightRanges{"rira", {2.42, 2.40}};
  SystematicDof<int> rebinFactors{"refa", {1, 2}};
  SystematicDof<int> bgFunctions{"bgfu", {2, 5}};

  TFile* fileMean = OpenFileWithNullptrCheck(meanCutVarFileName);
  TH1* histoMean = GetObjectWithNullptrCheck<TH1>(fileMean, "hCorrYieldsPrompt");

  const int nCtBins = histoMean->GetNbinsX();
  std::vector<TH1*> histoGets(nCtBins, nullptr);
  std::vector<TH1*> histoCounts(nCtBins, nullptr);
  for(int iCt=1; iCt<=nCtBins; ++iCt) {
    const double yield = histoMean->GetBinContent(iCt);
    for(auto* histos : {&histoGets, &histoCounts}) {
      const std::string histoName = (histos == &histoGets ? "yieldGet" : "yieldCount") + std::to_string(iCt);
      const double ctLo = histoMean->GetBinLowEdge(iCt);
      const double ctHi = histoMean->GetBinLowEdge(iCt+1);
      const std::string histoTitle = "Bin no." + std::to_string(iCt) + ", t#in(" + to_string_with_precision(ctLo, 1) + "#; " + to_string_with_precision(ctHi, 1) + ") ps";
      auto& histo = histos->at(iCt-1);
      histo = new TH1D(histoName.c_str(), histoTitle.c_str(), 100, 0.7*yield, 1.3*yield);
      histo->GetXaxis()->SetTitle("Corrected yield");
      histo->GetYaxis()->SetTitle("Entries");
    }
  }

  for(int iLr=0, nLrs=leftRanges.values_.size(); iLr<nLrs; ++iLr) {
    for(int iRr=0, nRrs=rightRanges.values_.size(); iRr<nRrs; ++iRr) {
      for(int iRf=0, nRfs=rebinFactors.values_.size(); iRf<nRfs; ++iRf) {
        for(int iBf=0, nBfs=bgFunctions.values_.size(); iBf<nBfs; ++iBf) {
          const std::string filePath = inputCommonPath + "/" +
                            leftRanges.prefix_ + "_" + to_string_with_precision(leftRanges.values_.at(iLr), 2) + "/" +
                            rightRanges.prefix_ + "_" + to_string_with_precision(rightRanges.values_.at(iRr), 2) + "/" +
                            rebinFactors.prefix_ + "_" + std::to_string(rebinFactors.values_.at(iRf)) + "/" +
                            bgFunctions.prefix_ + "_" + std::to_string(bgFunctions.values_.at(iBf));
          std::cout << "Info: Processing of " << filePath << "\n";
          auto ProcessYieldStrategy = [&] (const std::string& strategyName, std::vector<TH1*>& histoTargets) {
            TFile* fileYield{};
            TH1* histoYield{};
            try {
              fileYield = OpenFileWithNullptrCheck(filePath + "/" + strategyName + "/CutVarLc.merged.root");
              histoYield = GetObjectWithNullptrCheck<TH1>(fileYield, "hCorrYieldsPrompt");
            } catch(const std::exception&) {
              std::cout << "Info: Processing of " << filePath << "/" << strategyName << " is skipped due to missing file or histogram\n";
              return;
            }
            CheckHistogramsForXaxisIdentity(histoMean, histoYield);
            for(int iCt=1; iCt<=nCtBins; ++iCt) {
              const double yield = histoYield->GetBinContent(iCt);
              histoTargets.at(iCt-1)->Fill(yield);
            }
          };
          ProcessYieldStrategy("get", histoGets);
          ProcessYieldStrategy("count", histoCounts);
        } // bgFunctions
      } // rebinFactors
    } // rightRanges
  } // leftRanges

  TFile* fileOut = TFile::Open("corrected_yield_syst_qa.root", "recreate");
  for(int iCt=1; iCt<=nCtBins; ++iCt) {
    histoGets.at(iCt-1)->Write();
    histoCounts.at(iCt-1)->Write();
  }
  fileOut->Close();

  fileMean->Close();
}

int main(int argc, char* argv[]) {
  if (argc < 1) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./corrected_yield_syst_qa" << std::endl;
    exit(EXIT_FAILURE);
  }

  corrected_yield_syst_qa();

  return 0;
}
