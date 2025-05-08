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

std::vector<SimpleCut> PrepareBdtScoreCuts(const std::string& datatype, const std::string& direction, const std::vector<double>& vec, int precision=1);
std::vector<SimpleCut> PrepareBdtScoreCuts(const std::string& datatype, const std::string& direction, int nCuts, double loCut, double hiCut, int precision=1);

const std::vector<float> lifetimeRanges{0.2, 0.35, 0.5, 0.7, 0.9, 1.6};
const TAxis massAxis = {600, 1.98, 2.58};
const std::string massAxisTitle = "m_{pK#pi} (GeV/#it{c}^{2})";
const std::string massVarName = "hMass";

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
  man->SetVerbosityFrequency();
  man->Run(1000);
  man->Finish();
}

void MassQA(QA::Task& task) {
  auto TCuts = HelperFunctions::CreateRangeCuts(lifetimeRanges, "T_", branchName + ".KF_fT");

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
  auto TCuts = HelperFunctions::CreateRangeCuts(lifetimeRanges, "T_", branchName + ".fKFT", true);

  const std::string varNameInTree = branchName + ".fKFMassInv";

  const std::vector<SimpleCut> bdtBgCuts = PrepareBdtScoreCuts("bkg", "lt", 9, 0.1, 0.9);
  const std::vector<SimpleCut> bdtPromptCuts = PrepareBdtScoreCuts("prompt", "gt", 9, 0.1, 0.9);

  for(const auto& dt : datatypes) {
    for(const auto& bbc : bdtBgCuts) {
      for(const auto& bpc : bdtPromptCuts) {
        const std::string histoName = massVarName + "_" + bbc.GetTitle() + "_" + bpc.GetTitle();
        const QA::Axis histoQAAxis = {massAxisTitle, Variable::FromString(varNameInTree), massAxis};
        for(const auto& slc : TCuts) {
          std::string cutName = dt.GetTitle() + "/" + bbc.GetTitle();
          if(slc.GetTitle() != "alwaysTrue") cutName += "/" + slc.GetTitle();
          task.SetTopLevelDirName(cutName);
          Cuts* cutSlice = new Cuts(cutName, {dt, bbc, bpc, slc});
          task.AddH1(histoName, histoQAAxis, cutSlice);
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

std::vector<SimpleCut> PrepareBdtScoreCuts(const std::string& datatype, const std::string& direction, const std::vector<double>& vec, int precision) {
  if(datatype != "prompt" && datatype != "non_prompt" && datatype != "bkg") {
    throw std::runtime_error("mass_qa::PrepareBdtScoreCuts() - datatype must be either prompt, non_prompt or bkg");
  }
  if(direction != "gt" && direction != "lt") throw std::runtime_error("mass_qa::PrepareBdtScoreCuts() - direction must be either gt or lt");

  std::string dataLetter;
  if(datatype == "prompt") dataLetter = "P";
  if(datatype == "non_prompt") dataLetter = "NP";
  if(datatype == "bkg") dataLetter = "BG";

  std::vector<SimpleCut> result;
  for(const auto& value : vec) {
    const std::string sValue = HelperFunctions::ToStringWithPrecision(value, precision);
    const double lo = direction == "gt" ? value : -HugeNumber;
    const double hi = direction == "gt" ? HugeNumber : value;
    const std::string cutName = dataLetter + direction + sValue;
    result.emplace_back(RangeCut(branchName + "." + datatype + "_score", lo, hi, cutName));
  }
  return result;
}

std::vector<SimpleCut> PrepareBdtScoreCuts(const std::string& datatype, const std::string& direction, int nCuts, double loCut, double hiCut, int precision) {
  std::vector<double> vecCuts;
  const double step = (hiCut - loCut) / (nCuts - 1);
  for(int iCut=0; iCut<nCuts; iCut++) {
    const double value = loCut + iCut*step;
    vecCuts.emplace_back(value);
  }
  return PrepareBdtScoreCuts(datatype, direction, vecCuts, precision);
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
