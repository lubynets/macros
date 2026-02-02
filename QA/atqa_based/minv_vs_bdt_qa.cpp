//
// Created by oleksii on 23.01.2026.
//
#include "Task.hpp"

#include "AnalysisTree/HelperFunctions.hpp"
#include "AnalysisTree/TaskManager.hpp"

#include <iostream>
#include <string>
#include <vector>

using namespace AnalysisTree;

void MinvVsBdt(QA::Task& task);

void minv_vs_bdt_qa(const std::string& filelist, const int nEntries) {
  auto* man = TaskManager::GetInstance();

  auto* task = new QA::Task;
  task->SetOutputFileName("minv_vs_bdt_qa.root");

  MinvVsBdt(*task);

  man->AddTask(task);
  man->Init({filelist}, {"aTree"});
  man->SetVerbosityFrequency(100);
  man->Run(nEntries);
  man->Finish();
}

void MinvVsBdt(QA::Task& task) {
  const std::string sliceVarName{"PlainBranch.fLiteCt"};
  const std::vector<float> sliceVarRanges{0.006, 0.0105, 0.015, 0.021, 0.027, 0.048};
  const auto sliceCuts = HelperFunctions::CreateRangeCuts(sliceVarRanges, "Ct_", sliceVarName);

  for (const auto& sliceCut : sliceCuts) {
    const std::string& cutName = sliceCut.GetTitle();
    task.SetTopLevelDirName(cutName);
    Cuts* cut = new Cuts(cutName, {sliceCut});
    task.AddH2("hBdtVsMass", {"m_{pK#{pi}}(GeV/#it{c^2})", Variable::FromString("PlainBranch.fLiteM"), {600, 1.98, 2.58}}, {"BDT BG", Variable::FromString("PlainBranch.bkg_score"), {60, 0, 0.3}}, cut);
  }
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./minv_vs_bdt_qa filelistname (nEntries=ALL)\n";
    exit(EXIT_FAILURE);
  }

  const std::string filelistname = argv[1];
  const int nEntries = argc > 2 ? atoi(argv[2]) : -1;
  minv_vs_bdt_qa(filelistname, nEntries);

  return 0;
}