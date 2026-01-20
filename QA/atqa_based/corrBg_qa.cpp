//
// Created by oleksii on 23.12.2025.
//
#include "corrBg_qa.h"

#include "Task.hpp"

#include <AnalysisTree/HelperFunctions.hpp>
#include <AnalysisTree/TaskManager.hpp>

#include <iostream>
#include <string>

using namespace AnalysisTree;
using namespace HelperFunctions;

constexpr bool IsIncludeAllChannels = false;
constexpr bool IsApplyWeights = true;

const std::string recBranchName = "Candidates";
const std::pair<float, float> rapidityRanges{-0.8, 0.8};
std::vector<float> pTRanges = {1, 2, 3, 4, 5, 8, 12, 20};
const std::vector<float> bdtBgUpperValuesVsPt = {0.02, 0.02, 0.02, 0.02, 0.02, 0.04, 0.08};
const std::vector<float> lifetimeRanges = {0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.4, 3.6, 5.0}; const int lifetimeRangesPrecision = 1;

const TAxis massAxis = {600, 1.98, 2.58};
const std::string massAxisTitle = "m_{pK#pi} (GeV/#it{c}^{2})";

std::string GetCutName(size_t iCut, const std::vector<float>& ranges, const std::string& name, int precision);
std::string GetDecayFormula(const Decay& decay);

SimpleCut rapidityCut = RangeCut(Variable::FromString(recBranchName + ".fLiteY"), rapidityRanges.first, rapidityRanges.second);
Variable properLifetime("properLifetime", {{recBranchName, "fLiteCt"}}, [](const std::vector<double>& v) { return 100./2.997*v.at(0); });

enum KfSigBgStatus : int {
  kBackground = 0,
  kSignalPrompt = 1,
  kSignalNonPrompt = 2,
  kReflection = 3
};

enum : int {
  RunOnly = 0,
  RunAndMerge,
  MergeOnly,
  NModeRuns
};

void CorrBgQa(QA::Task& task) {
  if(bdtBgUpperValuesVsPt.size() != pTRanges.size() - 1) throw std::runtime_error("bdtUpperValuesVsPt.size() != pTRanges.size() - 1");

  for(int iPt=0, nPts=pTRanges.size()-1; iPt<nPts; ++iPt) {
    SimpleCut pTCut = RangeCut(recBranchName + ".fLitePt", pTRanges.at(iPt), pTRanges.at(iPt+1));
    SimpleCut bgBdtCut = RangeCut(recBranchName + ".fLiteMlScoreFirstClass", 0, bdtBgUpperValuesVsPt.at(iPt));
    const std::string pTCutName = "pT_" + HelperFunctions::ToStringWithPrecision(pTRanges.at(iPt), 0) + "_" + HelperFunctions::ToStringWithPrecision(pTRanges.at(iPt+1), 0);
    for(int iLifeTimeRange=0, nLifeTimeRanges=lifetimeRanges.size()-1; iLifeTimeRange<nLifeTimeRanges; ++iLifeTimeRange) {
      SimpleCut lifetimeCut = RangeCut(properLifetime, lifetimeRanges.at(iLifeTimeRange), lifetimeRanges.at(iLifeTimeRange+1));
      const std::string lifetimeCutName = "T_" + HelperFunctions::ToStringWithPrecision(lifetimeRanges.at(iLifeTimeRange), lifetimeRangesPrecision) + "_" + HelperFunctions::ToStringWithPrecision(lifetimeRanges.at(iLifeTimeRange+1), lifetimeRangesPrecision);
      task.SetTopLevelDirName(pTCutName + "/" + lifetimeCutName);

      for(const auto& decay : Decays) {
        if(!IsIncludeAllChannels && decay.is_bg_) continue;
        SimpleCut decayChannelCut = SimpleCut({recBranchName + ".fKFSigBgStatus", recBranchName + ".fLiteFlagMc"}, [&] (const std::vector<double>& var) {
          if(decay.is_bg_) {
            if(std::abs(var.at(1)) == decay.id_) {
              if(decay.id_ == LcToPKPi) return var.at(0) == kReflection;
              else return var.at(0) == kBackground;
            } else {
              return false;
            }
          } else {
            return var.at(0) == kSignalPrompt || var.at(0) == kSignalNonPrompt;
          }
        });
        const std::string decayFormula = GetDecayFormula(decay);
        Cuts* cutsDecay = new Cuts(decayFormula + "_" + pTCutName + "_" + lifetimeCutName, {rapidityCut, pTCut, bgBdtCut, lifetimeCut, decayChannelCut});
        Variable decayChannelWeight("decayChannelWeight", {{recBranchName, "ones"}}, [&] (std::vector<double>&) {
          if constexpr (IsApplyWeights) return decay.br_pdg_ / (decay.mother_.pythia_br_scaling_factor_ * decay.br_pythia_);
          else return 1.0;
        });
        task.AddH1("hMass_" + decayFormula, { massAxisTitle, Variable::FromString(recBranchName + ".fLiteM"), massAxis }, cutsDecay, decayChannelWeight);
      }

      SimpleCut bgChannelCut = SimpleCut({recBranchName + ".fKFSigBgStatus", recBranchName + ".fLiteFlagMc"}, [&] (const std::vector<double>& var) {
        for(const auto& decay : Decays) {
          if(std::abs(var.at(1)) == decay.id_) {
            if(decay.id_ == LcToPKPi) {
              if(var.at(0) == kReflection) return true;
            } else {
              if(var.at(0) == kBackground) return true;
            }
          }
        }
        return false;
      });

      Variable bgChannelWeight("bgChannelWeight", {{recBranchName, "fLiteFlagMc"}}, [&] (std::vector<double>& var) {
        if constexpr (IsApplyWeights) {
          for(const auto& decay : Decays) {
            if(!decay.is_bg_) continue;
            if(std::abs(decay.id_) == std::abs(var.at(0))) {
              return decay.br_pdg_ / (decay.mother_.pythia_br_scaling_factor_ * decay.br_pythia_);
            }
          }
          return 0.f;
        } else {
          return 1.f;
        }
      });

      Cuts* cutsBg = new Cuts("Bg_" + pTCutName + "_" + lifetimeCutName, {rapidityCut, pTCut, bgBdtCut, lifetimeCut, bgChannelCut});
      task.AddH1("hMass_bkgSum", { massAxisTitle, Variable::FromString(recBranchName + ".fLiteM"), massAxis }, cutsBg, bgChannelWeight);
    }
  }
}

void corrBg_qa(const std::string& fileInName, int modeRun) {
  const std::string& fileOutName = modeRun != MergeOnly ?  "corrBg_qa.root" : fileInName;

  if (modeRun != MergeOnly) {
    auto* man = TaskManager::GetInstance();
    auto* task = new QA::Task;
    task->SetOutputFileName(fileOutName);

    CorrBgQa(*task);

    man->AddTask(task);
    man->Init({fileInName}, {"aTree"});
    man->SetVerbosityFrequency(10);
    man->Run();
    man->Finish();
  }

  if (modeRun == RunOnly) return;

  const int nLowerPtBinsToExclude{0};
  pTRanges.erase(pTRanges.begin(), pTRanges.begin()+nLowerPtBinsToExclude);

  std::vector<std::string> pTCutNames, TCutNames;
  for(size_t iPt=0, nPts=pTRanges.size()-1; iPt<nPts; ++iPt) {
    pTCutNames.emplace_back(GetCutName(iPt, pTRanges, "pT", 0));
  }
  for(size_t iT=0, nTs=lifetimeRanges.size()-1; iT<nTs; ++iT) {
    TCutNames.emplace_back(GetCutName(iT, lifetimeRanges, "T", lifetimeRangesPrecision));
  }
  std::vector<Decay> DecaysWithBg{Decays};
  DecaysWithBg.emplace_back(Decay{-1, {}, "", -1.f, -1.f, true});

  TFile* fileOut = HelperFunctions::OpenFileWithNullptrCheck(fileOutName, "update");
  for(const auto& decay : DecaysWithBg) {
    if(!IsIncludeAllChannels && (decay.is_bg_ && decay.id_ != -1)) continue;
    for (const auto& tCut : TCutNames) {
      std::vector<std::string> histoNames;
      for (const auto& ptcn : pTCutNames) {
        histoNames.emplace_back(ptcn + "/" + tCut + "/hMass_" + GetDecayFormula(decay));
      }
    TH1* histoMerged = HelperFunctions::MergeHistograms(fileOut, histoNames);
    HelperFunctions::CD(fileOut, GetCutName(pTRanges.size()-1, pTRanges, "pT", 0) + "/" + tCut);
    histoMerged->Write(("hMass_" + GetDecayFormula(decay)).c_str());
    } // TCutNames
  } // DecaysWithBg
  fileOut->Close();
}

std::string GetCutName(const size_t iCut, const std::vector<float>& ranges, const std::string& name, const int precision) {
  std::pair<size_t, size_t> iPTMinMax = (iCut == ranges.size()-1) ? std::pair<size_t, size_t>{0, ranges.size()-1} : std::pair<size_t, size_t>{iCut, iCut+1};
  return name + "_" + HelperFunctions::ToStringWithPrecision(ranges.at(iPTMinMax.first), precision) + "_" + HelperFunctions::ToStringWithPrecision(ranges.at(iPTMinMax.second), precision);
}

std::string GetDecayFormula(const Decay& decay) {
  if(decay.id_ == -1) return "bkgSum";
  std::string decayFormula = decay.mother_.name_ + "To" + decay.daughters_;
  if(decay.id_ == LcToPKPi) decayFormula.append(decay.is_bg_ ? "Refl" : "Sig");

  return decayFormula;
}

int main(int argc, char* argv[]){
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./corrBg_qa fileInName (modeRun=RunOnly=0 [RunAndMerge=1, MergeOnly=2])" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileInName = argv[1];
  const int modeRun = argc > 2 ? std::stoi(argv[2]) : RunOnly;
  if(modeRun < 0 || modeRun >= NModeRuns) throw std::runtime_error("modeRun < 0 || modeRun >= NModeRuns");

  corrBg_qa(fileInName, modeRun);

  return 0;
}
