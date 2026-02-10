#include "Task.hpp"

#include <AnalysisTree/Constants.hpp>
#include <AnalysisTree/HelperFunctions.hpp>
#include <AnalysisTree/TaskManager.hpp>
#include <AnalysisTree/Variable.hpp>

using namespace AnalysisTree;

constexpr double MassLambdaCPlus{2.28646};
constexpr double LightSpeedCm2PS{0.0299792458};

void CtMcQa(QA::Task& task) {
  const std::array<SimpleCut, 2> promptnesses {
    EqualsCut(Variable::FromString("Candidates.fKFSigBgStatus"), 1, "Prompt"),
    EqualsCut(Variable::FromString("Candidates.fKFSigBgStatus"), 2, "NonPrompt")
  };

  const std::array<std::string, 2> mcStrategies{"Gen", "Kin"};

  SimpleCut ptCut = RangeCut(Variable::FromString("Candidates.fLitePt"), 3, 20);
  SimpleCut rapidityCut = RangeCut(Variable::FromString("Candidates.fLiteY"), -0.8, 0.8);
  constexpr std::array tRecRanges{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.4, 1.8, 2.4, 3.6, 5.0};

  Variable tRec("tRec", {{"Candidates", "fLiteCt"}}, []( std::vector<double>& v ) { return 100/2.997*v.at(0);});
  std::array<Variable, 2> tMc{
    Variable::FromString("Simulated.fSim_TDecay"),
    Variable("", {{"Simulated", "fSim_P"}, {"Simulated", "fSim_LDecay"}}, []( const std::vector<double>& v ) { return v.at(1) * MassLambdaCPlus / v.at(0) / LightSpeedCm2PS;})
  };
  std::array<Variable, 2> tDiff{
    Variable("tDiff", {{"Candidates", "fLiteCt"}, {"Simulated", "fSim_TDecay"}}, []( std::vector<double>& v ) { return 100/2.997*v.at(0) - v.at(1);}),
    Variable("tDiff", {{"Candidates", "fLiteCt"}, {"Simulated", "fSim_P"}, {"Simulated", "fSim_LDecay"}}, []( std::vector<double>& v ) { return 100/2.997*v.at(0) - (v.at(2) * MassLambdaCPlus / v.at(1) / LightSpeedCm2PS);})
  };

  const TAxis tAxis{400, 0, 2};
  const TAxis tDiffAxis{500, -0.5, 0.5};

  for(const auto& promptness : promptnesses) {
    Cuts* cutsGeneral = new Cuts("cutsGeneral", {ptCut, rapidityCut, promptness});
    for(int iMcStrategy=0; iMcStrategy<2; ++iMcStrategy) {
      task.SetTopLevelDirName(promptness.GetTitle() + "/" + mcStrategies.at(iMcStrategy));
      task.AddH2("timeCorr", {"t^{mc} (ps)", tMc.at(iMcStrategy), tAxis}, {"t^{rec} (ps)", tRec, tAxis}, cutsGeneral);
      task.AddH2("timeDiffVsMc", {"t^{mc} (ps)", tMc.at(iMcStrategy), tAxis}, {"t^{rec} - t^{mc} (ps)", tDiff.at(iMcStrategy), tDiffAxis}, cutsGeneral);
      task.AddH2("timeDiffVsRec", {"t^{rec} (ps)", tRec, tAxis}, {"t^{rec} - t^{mc} (ps)", tDiff.at(iMcStrategy), tDiffAxis}, cutsGeneral);
    }

    for(int iT=0, nTs = tRecRanges.size()-1; iT<nTs; ++iT) {
      SimpleCut tRceCut = RangeCut(tRec, tRecRanges.at(iT), tRecRanges.at(iT+1));
      Cuts* cutsTSilce = new Cuts("cutsTSilce", {ptCut, rapidityCut, tRceCut, promptness});
      const std::string hName = "timeDiff_T_" + HelperFunctions::ToStringWithPrecision(tRecRanges.at(iT), 1) + "_" + HelperFunctions::ToStringWithPrecision(tRecRanges.at(iT+1), 1);
      for(int iMcStrategy=0; iMcStrategy<2; ++iMcStrategy) {
        task.SetTopLevelDirName(promptness.GetTitle() + "/" + mcStrategies.at(iMcStrategy));
        task.AddH1(hName, {"t^{rec} - t^{mc} (ps)", tDiff.at(iMcStrategy), tDiffAxis}, cutsTSilce);
      }
    }
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
