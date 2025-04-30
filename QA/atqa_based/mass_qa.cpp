//
// Created by oleksii on 18.03.25.
//

#include "Task.hpp"

#include "AnalysisTree/HelperFunctions.hpp"
#include "AnalysisTree/TaskManager.hpp"

using namespace AnalysisTree;

void MassQA(QA::Task& task);
void MassQABdt(QA::Task& task);

const std::string branchName = "PlainBranch";

std::vector<SimpleCut> PrepareDataTypes(const std::string& varName);
std::vector<SimpleCut> datatypes = PrepareDataTypes(branchName + ".fKFSigBgStatus");

const std::vector<float> lifetimeRanges{0.2, 0.35, 0.5, 0.7, 0.9, 1.6};
const TAxis massAxis = {600, 1.98, 2.58};
const std::string massAxisTitle = "m_{pK#pi} (GeV/#it{c}^{2})";
const std::string massVarName = "Mass";

void mass_qa(const std::string& filelistname, bool isMc) {
  if(isMc) datatypes.pop_back();
  else     datatypes.erase(datatypes.begin(), datatypes.end()-1);

  auto* man = TaskManager::GetInstance();

  auto* task = new QA::Task;
  task->SetOutputFileName("mass_qa.root");

//  MassQA(*task);
  MassQABdt(*task);

  man->AddTask(task);
  man->Init({filelistname}, {"aTree"});
  man->SetVerbosityFrequency(100);
  man->Run(100000);
  man->Finish();
}

void MassQA(QA::Task& task) {
  auto TCuts = HelperFunctions::CreateRangeCuts(lifetimeRanges, "T_", branchName + ".KF_fT", 2);

  const std::string varNameInTree = branchName + ".KF_fMassInv";

  for(const auto& dt : datatypes) {
    task.SetTopLevelDirName(dt.GetTitle());
    Cuts* cutsTotal = new Cuts(dt.GetTitle() + "_total", {dt});
    task.AddH1(massVarName + "_" + dt.GetTitle(), {massAxisTitle, Variable::FromString(varNameInTree), massAxis}, cutsTotal);

    for(const auto& slc : TCuts) {
      Cuts* cutSlice = new Cuts(dt.GetTitle() + "_" + slc.GetTitle(), {dt, slc});
      task.AddH1(massVarName + "_" + dt.GetTitle() + "_" + slc.GetTitle(), {massAxisTitle, Variable::FromString(varNameInTree), massAxis}, cutSlice);
    } // TCuts
  } // datatypes
} // void MassQA()

void MassQABdt(QA::Task& task) {
  auto TCuts = HelperFunctions::CreateRangeCuts(lifetimeRanges, "T_", branchName + ".fKFT", 2);

  const std::string varNameInTree = branchName + ".fKFMassInv";

  std::vector<SimpleCut> bdtBgCuts;
  std::vector<SimpleCut> bdtPromptCuts;
  for(const auto& iv : std::vector<int>{2, 5, 8}) {
    bdtBgCuts.emplace_back(RangeCut(branchName + ".bkg_score", -HugeNumber, 1.*iv/10, "b" + std::to_string(iv)));
    bdtPromptCuts.emplace_back(RangeCut(branchName + ".prompt_score", 1.*iv/10, HugeNumber, "p" + std::to_string(iv)));
  }

  for(const auto& dt : datatypes) {
    for(const auto& bbc : bdtBgCuts) {
      for(const auto& bpc : bdtPromptCuts) {
        task.SetTopLevelDirName(dt.GetTitle());
        Cuts* cutsTotal = new Cuts(dt.GetTitle() + "_" + bbc.GetTitle() + bpc.GetTitle() + "_total", {dt, bbc, bpc});
        task.AddH1(massVarName + "_" + bbc.GetTitle() + bpc.GetTitle(), {massAxisTitle, Variable::FromString(varNameInTree), massAxis}, cutsTotal);

        for(const auto& slc : TCuts) {
          task.SetTopLevelDirName(dt.GetTitle() + "/" + slc.GetTitle());
          Cuts* cutSlice = new Cuts(dt.GetTitle() + "_" + bbc.GetTitle() + bpc.GetTitle() + "_" + slc.GetTitle(), {dt, bbc, bpc, slc});
          task.AddH1(massVarName + "_" + bbc.GetTitle() + bpc.GetTitle(), {massAxisTitle, Variable::FromString(varNameInTree), massAxis}, cutSlice);
        } // TCuts
      } // bdtPromptCuts
    } // bdtBgCuts
  } // datatypes
} // void MassQABdt

std::vector<SimpleCut> PrepareDataTypes(const std::string& varName) {
  std::vector<SimpleCut> result {
    RangeCut(varName, -0.1, 3.1, "all"),
    SimpleCut({varName}, [](std::vector<double> par){ return par[0] == 0 || par[0] == 1 || par[0] == 3; }, "all_wononprompt"),
    SimpleCut({varName}, [](std::vector<double> par){ return par[0] == 0 || par[0] == 2 || par[0] == 3; }, "all_woprompt"),
    SimpleCut({varName}, [](std::vector<double> par){ return par[0] == 0 || par[0] == 3; }, "background"),
    EqualsCut(varName, 1,"prompt"),
    EqualsCut(varName, 2,"nonprompt"),
    RangeCut(varName,  0.9, 2.1, "signal"),
    EqualsCut(varName, 3,"wrongswap"),
    EqualsCut(varName, -999, "data"),
  };

  return result;
}

int main(int argc, char* argv[]){
  if (argc < 3) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./mass_qa filelistname mcOrData" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string filelistname = argv[1];
  const bool isMc = strcmp(argv[2], "mc") == 0 ? true : strcmp(argv[2], "data") == 0 ? false : throw std::runtime_error("mass_qa::main(): argv[2] must be either 'mc' or 'data'");

  mass_qa(filelistname, isMc);

  return 0;
}
