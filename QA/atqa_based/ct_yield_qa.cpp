//
// Created by oleksii on 19.02.2026.
//
#include "Task.hpp"

#include <AnalysisTree/Constants.hpp>
#include <AnalysisTree/HelperFunctions.hpp>
#include <AnalysisTree/TaskManager.hpp>
#include <AnalysisTree/Variable.hpp>

#include <cstdlib>
#include <iostream>
#include <string>

using namespace AnalysisTree;

constexpr double MassLambdaCPlus{2.28646};
constexpr double LightSpeedCm2PS{0.0299792458};

const TAxis tAxis = {400, 0, 2};
const std::pair<double, double> rapidityRanges{-0.8, 0.8};
const std::pair<double, double> pTRanges{3., 20.};

struct Promptness {
  std::string name_;
  SimpleCut cut_rec_;
  SimpleCut cut_gen_;
};

const std::array<Promptness, 2> promptnesses {{
  {"prompt", EqualsCut("Candidates.fKFSigBgStatus", 1), EqualsCut("Generated.fGen_OriginMcGen", 1)},
  {"nonprompt", EqualsCut("Candidates.fKFSigBgStatus", 2), EqualsCut("Generated.fGen_OriginMcGen", 2)}
}};

SimpleCut ptCutCand = RangeCut(Variable::FromString("Candidates.fLitePt"), pTRanges.first, pTRanges.second);
SimpleCut rapidityCutCand = RangeCut(Variable::FromString("Candidates.fLiteY"), rapidityRanges.first, rapidityRanges.second);
SimpleCut ptCutGen = RangeCut(Variable::FromString("Generated.fGen_Pt"), pTRanges.first, pTRanges.second);
SimpleCut rapidityCutGen = RangeCut(Variable::FromString("Generated.fGen_Y"), rapidityRanges.first, rapidityRanges.second);

void CtYieldQa(QA::Task& task) {
  Variable tGen("tGen",
                {{"Generated", "fGen_P"}, {"Generated", "fGen_LDecay"}},
                [] (const std::vector<double>& v) {
                  return v.at(1) * MassLambdaCPlus / v.at(0) / LightSpeedCm2PS;
                });

  Variable tCand("tRec", {{"Candidates", "fLiteCt"}}, []( std::vector<double>& v ) { return v.at(0) / LightSpeedCm2PS;});

  Variable tSim("tGen",
                {{"Simulated", "fSim_P"}, {"Simulated", "fSim_LDecay"}},
                [] (const std::vector<double>& v) {
                  return v.at(1) * MassLambdaCPlus / v.at(0) / LightSpeedCm2PS;
                });

  for(const auto& promptness : promptnesses) {
    task.SetTopLevelDirName(promptness.name_.c_str());
    Cuts* cutsGen = new Cuts("cutsGen", {promptness.cut_gen_, ptCutGen, rapidityCutGen});
    task.AddH1("hGen", {"t^{gen} (ps)", tGen, tAxis}, cutsGen);

    Cuts* cutsRec = new Cuts("cutsRec", {promptness.cut_rec_, ptCutCand, rapidityCutCand});
    task.AddH1("hCand", {"t^{cand} (ps)", tCand, tAxis}, cutsRec);
    task.AddH1("hSim", {"t^{sim} (ps)", tSim, tAxis}, cutsRec);
  }
}

void ct_yield_qa(const std::string& filelist) {
  auto* man = TaskManager::GetInstance();

  auto* task = new QA::Task;
  task->SetOutputFileName("ct_yield_qa.root");

  CtYieldQa(*task);

  man->AddTask(task);
  man->Init({filelist}, {"aTree"});
  man->SetVerbosityFrequency(100);
  man->Run();
  man->Finish();
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./ct_yield_qa filelistname" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string filelistname = argv[1];
  ct_yield_qa(filelistname);

  return 0;
}
