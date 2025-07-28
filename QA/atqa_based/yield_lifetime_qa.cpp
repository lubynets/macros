//
// Created by oleksii on 02.05.25.
//

#include "Task.hpp"

#include "AnalysisTree/Constants.hpp"
#include "AnalysisTree/HelperFunctions.hpp"
#include "AnalysisTree/TaskManager.hpp"
#include "AnalysisTree/Variable.hpp"

#include <iostream>

using namespace AnalysisTree;

const std::vector<float> lifetimeRanges = {0.2, 0.35, 0.5, 0.7, 0.9, 1.6};
const std::vector<float> pTRanges = {0, 2, 5, 8, 12, 20};
const std::vector<float> bdtBgUpperValuesVsPt = {0.02, 0.02, 0.02, 0.05, 0.08};
const TAxis lifetimeAxis = {400, 0, 2};
const std::string lifetimeAxisTitle = "T (ps)";
const std::string genBranchName = "Generated";
const std::pair<float, float> rapidityRanges{-0.8, 0.8};

std::string GetPtCutName(size_t iPt) {
  if(iPt == pTRanges.size()-1) return "pT_" + HelperFunctions::ToStringWithPrecision(pTRanges.front(), 0) + "_" + HelperFunctions::ToStringWithPrecision(pTRanges.back(), 0);
  return "pT_" + HelperFunctions::ToStringWithPrecision(pTRanges.at(iPt), 0) + "_" + HelperFunctions::ToStringWithPrecision(pTRanges.at(iPt+1), 0);
}

enum BdtClass : short {
  kBackground = 0,
  kPrompt,
  kNonPrompt,
  nBdtClasses
};

const std::string recBranchName = "Candidates";
const std::vector<std::string> bdtClasses{"fLiteMlScoreFirstClass", "fLiteMlScoreSecondClass", "fLiteMlScoreThirdClass"};

// const std::string recBranchName = "PlainBranch";
//const std::vector<std::string> bdtClasses{"bkg_score", "prompt_score", "non_prompt_score"};

struct Promptness {
  std::string name_;
  SimpleCut cut_rec_;
  SimpleCut cut_gen_;
};

const std::vector<Promptness> promptnesses{
  {"prompt", EqualsCut(recBranchName + ".fKFSigBgStatus", 1), EqualsCut(genBranchName + ".fGen_OriginMcGen", 1)},
  {"nonprompt", EqualsCut(recBranchName + ".fKFSigBgStatus", 2), EqualsCut(genBranchName + ".fGen_OriginMcGen", 2)},
};

const std::vector<short> signalTargetsAtBdt {
//  kPrompt,
  kNonPrompt,
};

void FillYieldRec(QA::Task& task) {
  SimpleCut rapidityCut = RangeCut(Variable::FromString(recBranchName + ".fLiteY"), rapidityRanges.first, rapidityRanges.second);
  if (bdtBgUpperValuesVsPt.size() != pTRanges.size() - 1) throw std::runtime_error("bdtUpperValuesVsPt.size() != pTRanges.size() - 1");
  for(size_t iPt=0, nPts=pTRanges.size()-1; iPt<nPts; ++iPt) {
    const std::string pTCutName = GetPtCutName(iPt);
    SimpleCut bdtBgScoreCut = SimpleCut({recBranchName + ".fKFPt", recBranchName + "." + bdtClasses.at(kBackground)}, [=](const std::vector<double>& par) {
      return par[0] >= pTRanges.at(iPt) && par[0] < pTRanges.at(iPt + 1) && par[1] >= 0 && par[1] < bdtBgUpperValuesVsPt.at(iPt);
    }, pTCutName);
    for (const auto& kSignal : signalTargetsAtBdt) {
      const std::string signalShortcut = kSignal == kPrompt ? "P" : "NP";
      std::vector<float> bdtSignalLowerValues;
      for (int iB = 0; iB <= 99; iB++) {
        bdtSignalLowerValues.emplace_back(0.01 * iB);
      }

      std::vector<SimpleCut> bdtSigCuts;
      for (const auto& bslv : bdtSignalLowerValues) {
        bdtSigCuts.emplace_back(RangeCut(recBranchName + "." + bdtClasses.at(kSignal), bslv, 1, signalShortcut + "gt" + HelperFunctions::ToStringWithPrecision(bslv, 2)));
      }

      for (const auto& promptness : promptnesses) {
        task.SetTopLevelDirName("rec/" + promptness.name_ + "/" + pTCutName);
        for (const auto& bsc : bdtSigCuts) {
          Cuts* cutsRec = new Cuts(promptness.name_, {promptness.cut_rec_, bdtBgScoreCut, bsc, rapidityCut});
          const std::string histoName = "hT_" + bsc.GetTitle();
          task.AddH1(histoName, {lifetimeAxisTitle, Variable::FromString(recBranchName + ".fKFT"), lifetimeAxis}, cutsRec);
        }// bdtSigCuts
      }// promptnesses
    }// signalTargetsAtBdt
  } // pTs
} // FillYieldRec(QA::Task& task)

void FillYieldGen(QA::Task& task) {
  SimpleCut rapidityCut = RangeCut(Variable::FromString(genBranchName + ".fGen_Y"), rapidityRanges.first, rapidityRanges.second);
  for(size_t iPt=0, nPts=pTRanges.size()-1; iPt<nPts; ++iPt) {
    const std::string pTCutName = GetPtCutName(iPt);
    SimpleCut pTCut = RangeCut({genBranchName + ".fGen_Pt"}, pTRanges.at(iPt), pTRanges.at(iPt+1), pTCutName);
    for (const auto& promptness : promptnesses) {
      task.SetTopLevelDirName("gen/" + promptness.name_ + "/" + pTCutName);
      Cuts* cutsGen = new Cuts(promptness.name_ + "_" + pTCutName, {promptness.cut_gen_, pTCut});
      const std::string histoName = "hT";
      task.AddH1(histoName, {lifetimeAxisTitle, Variable::FromString(genBranchName + ".fGen_TDecay"), lifetimeAxis}, cutsGen);
    }// promptnesses
  } // pTs
} // FillYieldGen(QA::Task& task)

void yield_lifetime_qa(const std::string& filelistrec, const std::string& filelistgen) {
  auto* man = TaskManager::GetInstance();
  man->SetVerbosityFrequency(100);
  const std::string fileOutName = "yield_lifetime_qa.root";

  auto* taskRec = new QA::Task;
  taskRec->SetOutputFileName(fileOutName);

  FillYieldRec(*taskRec);

  man->AddTask(taskRec);
  man->Init({filelistrec}, {"aTree"});
  man->Run();
  man->Finish();

  man->ClearTasks();

  auto* taskGen = new QA::Task;
  taskGen->SetOutputFileName(fileOutName, "update");

  FillYieldGen(*taskGen);

  man->AddTask(taskGen);
  man->Init({filelistgen}, {"aTree"});
  man->Run();
  man->Finish();

  std::vector<std::string> pTCutNames;
  for(size_t iPt=0, nPts=pTRanges.size()-1; iPt<nPts; ++iPt) {
    pTCutNames.emplace_back(GetPtCutName(iPt));
  }
  TFile* fileOut = HelperFunctions::OpenFileWithNullptrCheck(fileOutName, "update");
  for(const auto& promptness : promptnesses) {
    std::vector<std::string> histoGenNames;
    for(const auto& ptcn : pTCutNames) {
      histoGenNames.emplace_back("gen/" + promptness.name_ + "/" + ptcn + "/hT");
    }
    TH1* histoGenMerged = HelperFunctions::MergeHistograms(fileOut, histoGenNames);
    HelperFunctions::CD(fileOut, "gen/" + promptness.name_ + "/" + GetPtCutName(pTRanges.size()-1));
    histoGenMerged->Write("hT");

    std::vector<float> bdtSignalLowerValues;
    for (int iB = 0; iB <= 99; iB++) {
      bdtSignalLowerValues.emplace_back(0.01 * iB);
    }
    for (const auto& kSignal : signalTargetsAtBdt) {
      const std::string signalShortcut = kSignal == kPrompt ? "P" : "NP";
      for (const auto& bslv : bdtSignalLowerValues) {
        std::vector<std::string> histoRecNames;
        for(const auto& ptcn : pTCutNames) {
          histoRecNames.emplace_back("rec/" + promptness.name_ + "/" + ptcn + "/hT_" + signalShortcut + "gt" + HelperFunctions::ToStringWithPrecision(bslv, 2));
        }
        TH1* histoRecMerged = HelperFunctions::MergeHistograms(fileOut, histoRecNames);
        HelperFunctions::CD(fileOut, "rec/" + promptness.name_ + "/" + GetPtCutName(pTRanges.size()-1));
        histoRecMerged->Write(("hT_" + signalShortcut + "gt" + HelperFunctions::ToStringWithPrecision(bslv, 2)).c_str());
      } // bdtSignalLowerValues
    } // signalTargetsAtBdt
  } // promptnesses
  fileOut->Close();
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./yield_lifetime_qa filelistrec (filelistgen=filelistrec)" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string filelistrec = argv[1];
  const std::string filelistgen = argc > 2 ? argv[2] : argv[1];
  yield_lifetime_qa(filelistrec, filelistgen);

  return 0;
}
