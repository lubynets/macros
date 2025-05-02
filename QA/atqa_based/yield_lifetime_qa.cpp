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
const TAxis lifetimeAxis = {400, 0, 2};
const std::string lifetimeAxisTitle = "T (ps)";

struct Promptness {
  std::string name_;
  SimpleCut cut_rec_;
  SimpleCut cut_gen_;
};

const std::vector<Promptness> promptnesses{
  {"prompt", EqualsCut("PlainBranch.fKFSigBgStatus", 1), EqualsCut("Generated.fGen_OriginMcGen", 1)},
  {"nonprompt", EqualsCut("PlainBranch.fKFSigBgStatus", 2), EqualsCut("Generated.fGen_OriginMcGen", 2)},
};

void FillYieldRec(QA::Task& task) {
  const std::string branchName = "PlainBranch";
  std::vector<SimpleCut> bdtBgCuts;
  std::vector<SimpleCut> bdtNonPromptCuts;
  const std::vector<float> cutValuesBg{0.003, 0.005, 0.008, 0.01, 0.012, 0.015, 0.02};
  for(const auto& cvb : cutValuesBg) {
    const std::string sCutValueBg = HelperFunctions::ToStringWithPrecision(cvb, 3);
    bdtBgCuts.emplace_back(RangeCut(branchName + ".bkg_score", -HugeNumber, cvb, "bg" + sCutValueBg));
  }
  for(int iv=1; iv<=10; iv++) {
    const double cutValueNp = 1.*iv/10;
    const std::string sCutValueNp = HelperFunctions::ToStringWithPrecision(cutValueNp, 1);
    bdtNonPromptCuts.emplace_back(RangeCut(branchName + ".non_prompt_score", -HugeNumber, cutValueNp, "np" + sCutValueNp));
  }

  for(const auto& promptness : promptnesses) {
    task.SetTopLevelDirName("rec/" + promptness.name_);
    for(const auto& bbc : bdtBgCuts) {
      for (const auto& bnpc : bdtNonPromptCuts) {
        Cuts* cutsRec = new Cuts(promptness.name_, {promptness.cut_rec_, bbc, bnpc});
        const std::string histoName = "hT_" + bbc.GetTitle() + "_" + bnpc.GetTitle();
        task.AddH1(histoName, {lifetimeAxisTitle, Variable::FromString("PlainBranch.fKFT"), lifetimeAxis}, cutsRec);
      } // bdtNonPromptCuts
    } // bdtBgCuts
  } // promptnesses
}

void FillYieldGen(QA::Task& task) {
  SimpleCut rapidityCut = RangeCut(Variable::FromString("Generated.fGen_Y"), -0.8, 0.8);
  for(const auto& promptness : promptnesses) {
    task.SetTopLevelDirName("gen/" + promptness.name_);
    Cuts* cutsGen = new Cuts(promptness.name_, {promptness.cut_gen_});
    const std::string histoName = "hT";
    task.AddH1(histoName, {lifetimeAxisTitle, Variable::FromString("Generated.fGen_TDecay"), lifetimeAxis}, cutsGen);
  } // promptnesses
}

void yield_lifetime_qa(const std::string& filelistrec, const std::string& filelistgen) {
  auto* man = TaskManager::GetInstance();
  man->SetVerbosityFrequency(100);

  auto* taskRec = new QA::Task;
  taskRec->SetOutputFileName("yield_lifetime_qa.root");

  FillYieldRec(*taskRec);

  man->AddTask(taskRec);
  man->Init({filelistrec}, {"aTree"});
  man->Run();
  man->Finish();

  man->ClearTasks();

  auto* taskGen = new QA::Task;
  taskGen->SetOutputFileName("yield_lifetime_qa.root", "update");

  FillYieldGen(*taskGen);

  man->AddTask(taskGen);
  man->Init({filelistgen}, {"aTree"});
  man->Run();
  man->Finish();
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./yield_lifetime_qa filelistrec filelistgen" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string filelistrec = argv[1];
  const std::string filelistgen = argv[2];
  yield_lifetime_qa(filelistrec, filelistgen);

  return 0;
}