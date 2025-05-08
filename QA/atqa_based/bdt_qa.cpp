//
// Created by oleksii on 28.04.25.
//
#include "Task.hpp"

#include "AnalysisTree/HelperFunctions.hpp"
#include "AnalysisTree/TaskManager.hpp"

#include <iostream>
#include <string>

using namespace AnalysisTree;

void BDTQA(QA::Task& task, const std::string& mcOrData);

void bdt_qa(const std::string& filelist, const std::string& mcOrData, int nEntries) {
  auto* man = TaskManager::GetInstance();

  auto* task = new QA::Task;
  task->SetOutputFileName("bdt_qa.root");

  BDTQA(*task, mcOrData);

  man->AddTask(task);
  man->Init({filelist}, {"aTree"});
  man->SetVerbosityFrequency(100);
  man->Run(nEntries);
  man->Finish();
}

void BDTQA(QA::Task& task, const std::string& mcOrData) {
  std::vector<SimpleCut> dataTypes;
  if(mcOrData == "mc") {
    dataTypes.emplace_back(EqualsCut("PlainBranch.fKFSigBgStatus", 1, "prompt"));
    dataTypes.emplace_back(EqualsCut("PlainBranch.fKFSigBgStatus", 2, "nonPrompt"));
  } else if (mcOrData == "data") {
    dataTypes.emplace_back(EqualsCut("PlainBranch.fKFSigBgStatus", -999, "background"));
  }

  auto pTCuts = HelperFunctions::CreateRangeCuts({0.f, 2.f, 5.f, 8.f, 12.f, 20.f}, "pT_", "PlainBranch.fKFPt");
  auto TCuts = HelperFunctions::CreateRangeCuts({0.2, 0.35, 0.5, 0.7, 0.9, 1.6}, "T_", "PlainBranch.fKFT", true);
  auto sliceCuts = HelperFunctions::MergeVectors(pTCuts, TCuts);

  const std::array<double, 4> sideBands{2.12, 2.20, 2.38, 2.42};
  SimpleCut sideBandsCut({"PlainBranch.fKFMassInv"}, [=] (const std::vector<double>& par) { return (par[0]>sideBands.at(0) && par[0]<sideBands.at(1)) || (par[0]>sideBands.at(2) && par[0]<sideBands.at(3)); }, "sidebands");

  const std::string histoName = "hBdt";

  const std::string xTitle = "bdt_{BG}";
  const Variable xVar = Variable::FromString("PlainBranch.bkg_score");
  const TAxis xAxis = {102, -0.01, 1.01};

  const std::string yTitle = "bdt_{Prompt}";
  const Variable yVar = Variable::FromString("PlainBranch.prompt_score");
  const TAxis yAxis = {102, -0.01, 1.01};

  for(const auto& dt : dataTypes) {
    for(const auto& sliceCut : sliceCuts) {
      std::string cutName = dt.GetTitle();
      if(sliceCut.GetTitle() != "alwaysTrue") cutName += "/" + sliceCut.GetTitle();
      task.SetTopLevelDirName(cutName);
      Cuts* cuts = new Cuts(cutName, {dt, sliceCut});
      if(dt.GetTitle() == "background") cuts->AddCut(sideBandsCut);
      task.AddH1(histoName, {xTitle, xVar, xAxis}, cuts);
      task.AddH2(histoName + "2D",{xTitle, xVar, xAxis}, {yTitle, yVar, yAxis}, cuts);
    } // sliceCuts
  } // dataTypes
}

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./bdt_qa filelistname mcOrData (nEntries=ALL)" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string filelistname = argv[1];
  const std::string mcOrData = argv[2];
  if (mcOrData != "mc" && mcOrData != "data") throw std::runtime_error("bdt_qa::main(): mcOrData must be either 'mc' or 'data'");
  const int nEntries = argc > 3 ? atoi(argv[3]) : -1;
  bdt_qa(filelistname, mcOrData, nEntries);

  return 0;
}