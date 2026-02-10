//
// Created by oleksii on 10.02.2026.
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

void CtGenQa(QA::Task& task) {
  Variable properDecayTimeGenerator(Variable::FromString("Generated.fGen_TDecay"));
  Variable properDecayTimeKinematic("properDecayTimeKinematic",
                                    {{"Generated", "fGen_P"}, {"Generated", "fGen_LDecay"}},
                                    [] (const std::vector<double>& v) {
                                      return v.at(1) * MassLambdaCPlus / v.at(0) / LightSpeedCm2PS;
                                    });
  Variable properDecayTimeDiffKinGen("properDecayTimeDiffKinGen",
                                     {{"Generated", "fGen_P"}, {"Generated", "fGen_LDecay"}, {"Generated", "fGen_TDecay"}},
                                     [] (const std::vector<double>& v) {
                                       return v.at(1) * MassLambdaCPlus / v.at(0) / LightSpeedCm2PS - v.at(2);
                                     });

  const std::array<SimpleCut, 2> promptnesses {
    EqualsCut(Variable::FromString("Generated.fGen_OriginMcGen"),1, "Prompt"),
    EqualsCut(Variable::FromString("Generated.fGen_OriginMcGen"),2, "NonPrompt")
  };

  const TAxis tAxis{400, 0, 2};
  const TAxis tDiffAxis{2000, -0.5, 0.5};

  task.SetTopLevelDirName("");
  for(const auto& promptness : promptnesses) {
    Cuts* cuts = new Cuts("cutsPromptness", {promptness});
    task.AddH2("hCorr" + promptness.GetTitle(), {"t^{gen} (ps)", properDecayTimeGenerator, tAxis}, {"t^{kin} (ps)", properDecayTimeKinematic, tAxis}, cuts);
    task.AddH1("hDiff" + promptness.GetTitle(), {"t^{kin} - t^{gen} (ps)", properDecayTimeDiffKinGen, tDiffAxis}, cuts);
  }
}

void ct_gen_qa(const std::string& filelist) {
  auto* man = TaskManager::GetInstance();

  auto* task = new QA::Task;
  task->SetOutputFileName("ct_gen_qa.root");

  CtGenQa(*task);

  man->AddTask(task);
  man->Init({filelist}, {"aTree"});
  man->SetVerbosityFrequency(100);
  man->Run();
  man->Finish();
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./ct_gen_qa filelistname" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string filelistname = argv[1];
  ct_gen_qa(filelistname);

  return 0;
}