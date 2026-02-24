#include "Task.hpp"

#include <AnalysisTree/Constants.hpp>
#include <AnalysisTree/HelperFunctions.hpp>
#include <AnalysisTree/TaskManager.hpp>
#include <AnalysisTree/Variable.hpp>

#include <cmath>

using namespace AnalysisTree;

constexpr double MassLambdaCPlus{2.28646};
constexpr double LightSpeedCm2PS{0.0299792458};

void CtMcQa(QA::Task& task) {
  const std::array<SimpleCut, 2> promptnesses {
    EqualsCut(Variable::FromString("Candidates.fKFSigBgStatus"), 1, "Prompt"),
    EqualsCut(Variable::FromString("Candidates.fKFSigBgStatus"), 2, "NonPrompt")
  };

  const std::array<std::string, 2> signStrategies{"unsigned", "signed"};

  SimpleCut ptCut = RangeCut(Variable::FromString("Candidates.fLitePt"), 3, 20);
  SimpleCut rapidityCut = RangeCut(Variable::FromString("Candidates.fLiteY"), -0.8, 0.8);

  std::array<Variable, 2> tRec {
    Variable("tRec", {{"Candidates", "fLiteCt"}}, []( std::vector<double>& v ) { return v.at(0) / LightSpeedCm2PS;}),

    Variable("tRec", {
      {"Candidates", "fLitePosX"}, // 0
      {"Candidates", "fLitePosY"}, // 1
      {"Candidates", "fLitePosZ"}, // 2
      {"Candidates", "fKFX"},      // 3
      {"Candidates", "fKFY"},      // 4
      {"Candidates", "fKFZ"},      // 5
      {"Candidates", "fLitePt"},   // 6
      {"Candidates", "fLiteEta"},  // 7
      {"Candidates", "fLitePhi"}   // 8
    }, []( std::vector<double>& v ) {
      const auto dx = v[3] - v[0];
      const auto dy = v[4] - v[1];
      const auto dz = v[5] - v[2];
      const auto px = v[6] * std::cos(v[8]);
      const auto py = v[6] * std::sin(v[8]);
      const auto pz = v[6] * std::sinh(v[7]);
      const auto lp = dx*px + dy*py + dz*pz;
      const auto p2 = px*px + py*py + pz*pz;
      return MassLambdaCPlus * lp / LightSpeedCm2PS / p2;
    })
  };

  Variable tMc("", {{"Simulated", "fSim_P"}, {"Simulated", "fSim_LDecay"}}, []( const std::vector<double>& v ) { return v.at(1) * MassLambdaCPlus / v.at(0) / LightSpeedCm2PS;});

  std::array<Variable, 2> tDiff {
    Variable("tDiff", {{"Candidates", "fLiteCt"}, {"Simulated", "fSim_P"}, {"Simulated", "fSim_LDecay"}}, []( std::vector<double>& v ) { return v.at(0) / LightSpeedCm2PS - (v.at(2) * MassLambdaCPlus / v.at(1) / LightSpeedCm2PS);}),

    Variable("tDiff", {
      {"Candidates", "fLitePosX"}, // 0
      {"Candidates", "fLitePosY"}, // 1
      {"Candidates", "fLitePosZ"}, // 2
      {"Candidates", "fKFX"},      // 3
      {"Candidates", "fKFY"},      // 4
      {"Candidates", "fKFZ"},      // 5
      {"Candidates", "fLitePt"},   // 6
      {"Candidates", "fLiteEta"},  // 7
      {"Candidates", "fLitePhi"},  // 8
      {"Simulated", "fSim_P"},     // 9
      {"Simulated", "fSim_LDecay"} // 10
    }, []( std::vector<double>& v ) {
      const auto dx = v[3] - v[0];
      const auto dy = v[4] - v[1];
      const auto dz = v[5] - v[2];
      const auto px = v[6] * std::cos(v[8]);
      const auto py = v[6] * std::sin(v[8]);
      const auto pz = v[6] * std::sinh(v[7]);
      const auto lp = dx*px + dy*py + dz*pz;
      const auto p2 = px*px + py*py + pz*pz;
      return MassLambdaCPlus / LightSpeedCm2PS * (lp / p2 - v[10] / v[9]);
    })
  };

  const TAxis tAxis{600, -1, 2};
  const TAxis tDiffAxis{500, -0.5, 0.5};

  for(const auto& promptness : promptnesses) {
    Cuts* cutsGeneral = new Cuts("cutsGeneral", {ptCut, rapidityCut, promptness});
    for(int iSign=0, nSigns=signStrategies.size(); iSign<nSigns; ++iSign) {
      task.SetTopLevelDirName(promptness.GetTitle() + "/" + signStrategies.at(iSign));
      task.AddH2("timeCorr", {"t^{mc} (ps)",tMc, tAxis}, {"t^{rec} (ps)", tRec.at(iSign), tAxis}, cutsGeneral);
      task.AddH2("timeDiffVsMc", {"t^{mc} (ps)",tMc, tAxis}, {"t^{rec} - t^{mc} (ps)", tDiff.at(iSign), tDiffAxis}, cutsGeneral);
      task.AddH2("timeDiffVsRec", {"t^{rec} (ps)",tRec.at(iSign), tAxis}, {"t^{rec} - t^{mc} (ps)", tDiff.at(iSign), tDiffAxis}, cutsGeneral);
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
