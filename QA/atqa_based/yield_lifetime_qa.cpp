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
const TAxis lifetimeAxis = {200, 0, 2};
const std::string lifetimeAxisTitle = "T (ps)";
const std::string genBranchName = "Generated";
const std::pair<float, float> rapidityRanges{-0.8, 0.8};

std::string GetPtCutName(size_t iPt) {
  std::pair<size_t, size_t> iPTMinMax = (iPt == pTRanges.size()-1) ? std::pair<size_t, size_t>{0, pTRanges.size()-1} : std::pair<size_t, size_t>{iPt, iPt+1};
  return "pT_" + HelperFunctions::ToStringWithPrecision(pTRanges.at(iPTMinMax.first), 0) + "_" + HelperFunctions::ToStringWithPrecision(pTRanges.at(iPTMinMax.second), 0);
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

const std::vector<std::string> weightsPresences{"", "_W"};

const short kSignal = kNonPrompt; const std::string signalShortcut = "NP";

Variable properLifetime("properLifetime", {{recBranchName, "fLiteCt"}}, [](const std::vector<double>& v) { return 100./2.99792458*v.at(0); });

void FillYieldRec(QA::Task& task, const TH1* hWeights) {
  SimpleCut rapidityCut = RangeCut(Variable::FromString(recBranchName + ".fLiteY"), rapidityRanges.first, rapidityRanges.second);
  Variable ptWeigth("ptWeigth", {{recBranchName, "fKFPt"}}, [=] (std::vector<double>& par) {
    return HelperFunctions::InterpolateTH1SuppressWarning(hWeights, par[0]);
  } );
  if (bdtBgUpperValuesVsPt.size() != pTRanges.size() - 1) throw std::runtime_error("bdtUpperValuesVsPt.size() != pTRanges.size() - 1");
  for(size_t iPt=0, nPts=pTRanges.size()-1; iPt<nPts; ++iPt) {
    const std::string pTCutName = GetPtCutName(iPt);
    SimpleCut bdtBgScoreCut = SimpleCut({recBranchName + ".fKFPt", recBranchName + "." + bdtClasses.at(kBackground)}, [=](const std::vector<double>& par) {
      return par[0] >= pTRanges.at(iPt) && par[0] < pTRanges.at(iPt + 1) && par[1] >= 0 && par[1] < bdtBgUpperValuesVsPt.at(iPt);
    }, pTCutName);
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
        task.AddH1(histoName, {lifetimeAxisTitle, properLifetime, lifetimeAxis}, cutsRec);
        if(promptness.name_ == "prompt") task.AddH1(histoName + "_W", {lifetimeAxisTitle, properLifetime, lifetimeAxis}, cutsRec, ptWeigth);
      }// bdtSigCuts
    }// promptnesses
  } // pTs
} // FillYieldRec(QA::Task& task)

void FillYieldGen(QA::Task& task, const TH1* hWeights) {
  SimpleCut rapidityCut = RangeCut(Variable::FromString(genBranchName + ".fGen_Y"), rapidityRanges.first, rapidityRanges.second);
  Variable ptWeigth("ptWeigth", {{genBranchName, "fGen_Pt"}}, [=] (std::vector<double>& par) {
    return HelperFunctions::InterpolateTH1SuppressWarning(hWeights, par[0]);
  } );
  for(size_t iPt=0, nPts=pTRanges.size()-1; iPt<nPts; ++iPt) {
    const std::string pTCutName = GetPtCutName(iPt);
    SimpleCut pTCut = RangeCut({genBranchName + ".fGen_Pt"}, pTRanges.at(iPt), pTRanges.at(iPt+1), pTCutName);
    for (const auto& promptness : promptnesses) {
      task.SetTopLevelDirName("gen/" + promptness.name_ + "/" + pTCutName);
      Cuts* cutsGen = new Cuts(promptness.name_ + "_" + pTCutName, {promptness.cut_gen_, pTCut, rapidityCut});
      const std::string histoName = "hT";
      task.AddH1(histoName, {lifetimeAxisTitle, Variable::FromString(genBranchName + ".fGen_TDecay"), lifetimeAxis}, cutsGen);
      task.AddH1(histoName + "_W", {lifetimeAxisTitle, Variable::FromString(genBranchName + ".fGen_TDecay"), lifetimeAxis}, cutsGen, ptWeigth);
    }// promptnesses
  } // pTs
} // FillYieldGen(QA::Task& task)

void yield_lifetime_qa(const std::string& filelistrec, const std::string& filelistgen, const std::string& filePtWeights) {
  auto* man = TaskManager::GetInstance();
  man->SetVerbosityFrequency(100);
  const std::string fileOutName = "yield_lifetime_qa.root";

  TFile* fileWeight = HelperFunctions::OpenFileWithNullptrCheck(filePtWeights, "read");
  TH1* histoWeight = HelperFunctions::GetObjectWithNullptrCheck<TH1>(fileWeight, "histoWeight_pT_0_20");

  auto* taskRec = new QA::Task;
  taskRec->SetOutputFileName(fileOutName);

  FillYieldRec(*taskRec, histoWeight);

  man->AddTask(taskRec);
  man->Init({filelistrec}, {"aTree"});
  man->Run();
  man->Finish();

  man->ClearTasks();

  auto* taskGen = new QA::Task;
  taskGen->SetOutputFileName(fileOutName, "update");

  FillYieldGen(*taskGen, histoWeight);

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
    for(const auto& weightPresence : weightsPresences) {
      if(promptness.name_ == "nonprompt" && weightPresence == "_W") continue;
      std::vector<std::string> histoGenNames;
      for (const auto& ptcn : pTCutNames) {
        histoGenNames.emplace_back("gen/" + promptness.name_ + "/" + ptcn + "/hT" + weightPresence);
      }
      TH1* histoGenMerged = HelperFunctions::MergeHistograms(fileOut, histoGenNames);
      HelperFunctions::CD(fileOut, "gen/" + promptness.name_ + "/" + GetPtCutName(pTRanges.size() - 1));
      histoGenMerged->Write(("hT" + weightPresence).c_str());
    } // weightPresences

    std::vector<float> bdtSignalLowerValues;
    for (int iB = 0; iB <= 99; iB++) {
      bdtSignalLowerValues.emplace_back(0.01 * iB);
    }
    for (const auto& bslv : bdtSignalLowerValues) {
      for(const auto& weightPresence : weightsPresences) {
        if(promptness.name_ == "nonprompt" && weightPresence == "_W") continue;
        std::vector<std::string> histoRecNames;
        for (const auto& ptcn : pTCutNames) {
          histoRecNames.emplace_back("rec/" + promptness.name_ + "/" + ptcn + "/hT_" + signalShortcut + "gt" + HelperFunctions::ToStringWithPrecision(bslv, 2) + weightPresence);
        }
        TH1* histoRecMerged = HelperFunctions::MergeHistograms(fileOut, histoRecNames);
        HelperFunctions::CD(fileOut, "rec/" + promptness.name_ + "/" + GetPtCutName(pTRanges.size() - 1));
        histoRecMerged->Write(("hT_" + signalShortcut + "gt" + HelperFunctions::ToStringWithPrecision(bslv, 2) + weightPresence).c_str());
      } // weightPresences
    } // bdtSignalLowerValues
  } // promptnesses
  fileOut->Close();
  fileWeight->Close();
}

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./yield_lifetime_qa filelistrec filePtWeights" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string filelistrec = argv[1];
  const std::string filelistgen = argv[1];
  const std::string filePtWeights = argv[2];
  yield_lifetime_qa(filelistrec, filelistgen, filePtWeights);

  return 0;
}
