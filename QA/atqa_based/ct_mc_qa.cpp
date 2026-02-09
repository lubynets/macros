#include "Task.hpp"

#include <AnalysisTree/Constants.hpp>
#include <AnalysisTree/HelperFunctions.hpp>
#include <AnalysisTree/TaskManager.hpp>
#include <AnalysisTree/Variable.hpp>

using namespace AnalysisTree;

void CtMcQa(QA::Task& task) {
  SimpleCut ptCut = RangeCut(Variable::FromString("Candidates.fLitePt"), 3, 20);
  SimpleCut rapidityCut = RangeCut(Variable::FromString("Candidates.fLiteY"), -0.8, 0.8);
  constexpr std::array tRecRanges{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.4, 1.8, 2.4, 3.6, 5.0};

  Variable tRec("tRec", {{"Candidates", "fLiteCt"}}, []( std::vector<double>& v ) { return 100/2.997*v.at(0);});
  Variable tMc(Variable::FromString("Simulated.fSim_TDecay"));
  Variable tDiff("tDiff", {{"Candidates", "fLiteCt"}, {"Simulated", "fSim_TDecay"}}, []( std::vector<double>& v ) { return 100/2.997*v.at(0) - v.at(1);});

  const TAxis tAxis{400, 0, 2};
  const TAxis tDiffAxis{500, -0.5, 0.5};

  Cuts* cutsGeneral = new Cuts("cutsGeneral", {ptCut, rapidityCut});
  task.SetTopLevelDirName("");
  task.AddH2("timeCorr", {"t^{mc} (ps)", tMc, tAxis}, {"t^{rec} (ps)", tRec, tAxis}, cutsGeneral);
  task.AddH2("timeDiffVsMc", {"t^{mc} (ps)", tMc, tAxis}, {"t^{rec} - t^{mc} (ps)", tDiff, tDiffAxis}, cutsGeneral);
  task.AddH2("timeDiffVsRec", {"t^{rec} (ps)", tRec, tAxis}, {"t^{rec} - t^{mc} (ps)", tDiff, tDiffAxis}, cutsGeneral);

  for(int iT=0, nTs = tRecRanges.size()-1; iT<nTs; ++iT) {
    SimpleCut tRceCut = RangeCut(tRec, tRecRanges.at(iT), tRecRanges.at(iT+1));
    Cuts* cutsTSilce = new Cuts("cutsTSilce", {ptCut, rapidityCut, tRceCut});
    const std::string hName = "timeDiff_T_" + HelperFunctions::ToStringWithPrecision(tRecRanges.at(iT), 1) + "_" + HelperFunctions::ToStringWithPrecision(tRecRanges.at(iT+1), 1);
    task.AddH1(hName, {"t^{rec} - t^{mc} (ps)", tDiff, tDiffAxis}, cutsTSilce);
  }
}

void ct_mc_qa(const std::string& filelist) {
  auto* man = TaskManager::GetInstance();

  auto* task = new QA::Task;
  task->SetOutputFileName("ct_mc_qa.root");

  CtMcQa(*task);

  man->AddTask(task);
  man->Init({filelist}, {"aTree"});
  man->SetVerbosityFrequency(100);
  man->Run();
  man->Finish();
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./ct_mc_qa filelistname" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string filelistname = argv[1];
  ct_mc_qa(filelistname);

  return 0;
}
