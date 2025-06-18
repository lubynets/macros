//
// Created by oleksii on 18.03.25.
//

#include "Task.hpp"

#include "AnalysisTree/HelperFunctions.hpp"
#include "AnalysisTree/TaskManager.hpp"

using namespace AnalysisTree;

const std::vector<float> lifetimeRanges = {0.2, 0.35, 0.5, 0.7, 0.9, 1.6};
const std::vector<float> pTRanges = {0, 2, 5, 8, 12, 20};
const TAxis massAxis = {600, 1.98, 2.58};
const std::string massAxisTitle = "m_{pK#pi} (GeV/#it{c}^{2})";
const std::string massVarName = "hMass";
const std::pair<float, float> rapidityRanges{-0.8, 0.8};

enum BdtClass : short {
  kBackground = 0,
  kPrompt,
  kNonPrompt,
  nBdtClasses
};

const std::string recBranchName = "PlainBranch";
const std::vector<std::string> bdtClasses{"fLiteMlScoreFirstClass", "fLiteMlScoreSecondClass", "fLiteMlScoreThirdClass"};

auto TCuts = HelperFunctions::CreateRangeCuts(lifetimeRanges, "T_", recBranchName + ".fKFT");
SimpleCut rapidityCut = RangeCut(Variable::FromString(recBranchName + ".fLiteY"), rapidityRanges.first, rapidityRanges.second);

void MassQA(QA::Task& task);
void MassQABdt(QA::Task& task, const std::string& signalTargetAtBdt);

std::vector<SimpleCut> PrepareDataTypes(const std::string& varName);
std::vector<SimpleCut> datatypes = PrepareDataTypes(recBranchName + ".fKFSigBgStatus");

void mass_qa(const std::string& filelistname, bool isMc, const std::string& signalTargetAtBdt) {
  if(isMc) datatypes.pop_back();
  else     datatypes.erase(datatypes.begin(), datatypes.end()-1);

  auto* man = TaskManager::GetInstance();

  auto* task = new QA::Task;
  task->SetOutputFileName("mass_qa.root");

//  MassQA(*task);
  MassQABdt(*task, signalTargetAtBdt);

  man->AddTask(task);
  man->Init({filelistname}, {"aTree"});
  man->SetVerbosityFrequency(100);
  man->Run();
  man->Finish();
}

void MassQA(QA::Task& task) {
  const std::string varNameInTree = recBranchName + ".KF_fMassInv";

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

void MassQABdt(QA::Task& task, const std::string& signalTargetAtBdt) {
  const short kSignal = signalTargetAtBdt == "prompt" ? kPrompt : kNonPrompt;
  const std::string signalShortcut = signalTargetAtBdt == "prompt" ? "P" : "NP";

  const std::vector<float> bdtBgUpperValuesVsPt = {0.02, 0.02, 0.02, 0.05, 0.08};
  if(bdtBgUpperValuesVsPt.size() != pTRanges.size() - 1) throw std::runtime_error("bdtUpperValuesVsPt.size() != pTRanges.size() - 1");
  SimpleCut bdtBgScoreCut = SimpleCut({recBranchName + ".fKFPt", recBranchName + "." + bdtClasses.at(kBackground)}, [=] (const std::vector<double>& par) {
    bool ok{false};
    for(int iPt=0; iPt<bdtBgUpperValuesVsPt.size(); iPt++) {
      ok |= (par[0] >= pTRanges.at(iPt) && par[0] < pTRanges.at(iPt+1) && par[1] >= 0 && par[1] < bdtBgUpperValuesVsPt.at(iPt));
    }
    return ok;
  });

  std::vector<float> bdtSignalLowerValues;
  for(int iB=0; iB<=50; iB++) {
    bdtSignalLowerValues.emplace_back(0.02 * iB);
  }

  std::vector<SimpleCut> bdtSigCuts;
  for(const auto& bslv : bdtSignalLowerValues) {
    bdtSigCuts.emplace_back(RangeCut(recBranchName + "." + bdtClasses.at(kSignal), bslv, 1, signalShortcut + "gt" + HelperFunctions::ToStringWithPrecision(bslv, 2)));
  }

  for(const auto& dt : datatypes) {
    for(const auto& tCut : TCuts) {
      task.SetTopLevelDirName(dt.GetTitle() + "/" + tCut.GetTitle());
      for (const auto& bsc : bdtSigCuts) {
        Cuts* cutsRec = new Cuts(dt.GetTitle(), {bdtBgScoreCut, bsc, rapidityCut, tCut});
        const std::string histoName = "hM_" + bsc.GetTitle();
        task.AddH1(histoName, {massAxisTitle, Variable::FromString(recBranchName + ".fKFMassInv"), massAxis}, cutsRec);
      } // bdtNonPromptCuts
    } // TCuts
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
    std::cout << " ./mass_qa filelistname mcOrData (signalTargetAtBdt=prompt)" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string filelistname = argv[1];
  const bool isMc = strcmp(argv[2], "mc") == 0 ? true : strcmp(argv[2], "data") == 0 ? false : throw std::runtime_error("mass_qa::main(): argv[2] must be either 'mc' or 'data'");
  const std::string signalTargetAtBdt = argc > 3 ? argv[3] : "prompt";
  if(signalTargetAtBdt != "prompt" && signalTargetAtBdt != "nonprompt") {
    throw std::runtime_error("mass_qa::main(): signalTargetAtBdt must be either 'prompt' or 'nonprompt'");
  }

  mass_qa(filelistname, isMc, signalTargetAtBdt);

  return 0;
}
