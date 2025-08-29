#include "AnalysisTree/HelperFunctions.hpp"
#include "AnalysisTree/TaskManager.hpp"
#include "AnalysisTree/Variable.hpp"
#include "Task.hpp"

using namespace AnalysisTree;

void KfVsDca(QA::Task& task, bool isMc);

void kf_vs_dca_qa(const std::string& filelist, bool isMc){
  auto* man = TaskManager::GetInstance();

  auto* task = new QA::Task;
  task->SetOutputFileName("kf_vs_dca_qa.root");

  KfVsDca(*task, isMc);

  man->AddTask(task);
  man->Init({filelist}, {"aTree"});
  man->SetVerbosityFrequency(100);
  man->Run();
  man->Finish();
}

void KfVsDca(QA::Task& task, bool isMc) {
  Variable tDCA("tDCA", {{"Candidates", "fLiteCt"}}, []( std::vector<double>& v ) { return 100/2.997*v.at(0);});
  task.AddH2("timeCorr",{"t_{proper}^{DCAFitter} (fs)", tDCA, {500, -0.5, 2}}, {"t_{proper}^{KF} (fs)", Variable::FromString("Candidates.fKFT"), {500, -0.5, 2}});
  Variable tDiff("tDiff", {{"Candidates", "fLiteCt"}, {"Candidates", "fKFT"}}, []( std::vector<double>& v ) { return 100/2.997*v.at(0) - v.at(1); });
  task.AddH1("timeDiff", {"t_{proper}^{DCAFitter} - t_{proper}^{KF} (fs)", tDiff, {3000, -0.1, 0.8}});

  if(isMc) {
    task.AddH2("timeCorrDCAWithMc",{"t_{proper}^{DCAFitter} (fs)", tDCA, {500, -0.5, 2}}, {"t_{proper}^{MC} (fs)", Variable::FromString("Simulated.fSim_TDecay"), {500, -0.5, 2}});
    task.AddH2("timeCorrKFWithMc",{"t_{proper}^{KF} (fs)", Variable::FromString("Candidates.fKFT"), {500, -0.5, 2}}, {"t_{proper}^{MC} (fs)", Variable::FromString("Simulated.fSim_TDecay"), {500, -0.5, 2}});
    Variable tDiffDCAWithMc("tDiffDCAWithMc", {{"Candidates", "fLiteCt"}, {"Simulated", "fSim_TDecay"}}, []( std::vector<double>& v ) { return 100/2.997*v.at(0) - v.at(1); });
    Variable tDiffKFWithMc("tDiffKFWithMc", {{"Candidates", "fKFT"}, {"Simulated", "fSim_TDecay"}}, []( std::vector<double>& v ) { return v.at(0) - v.at(1); });
    task.AddH1("timeDiffDCAWithMc", {"t_{proper}^{DCAFitter} - t_{proper}^{MC} (fs)", tDiffDCAWithMc, {3000, -0.9, 0.9}});
    task.AddH1("timeDiffKFWithMc", {"t_{proper}^{KF} - t_{proper}^{MC} (fs)", tDiffKFWithMc, {3000, -0.9, 0.9}});
  }

  task.AddH2("massCorr",{"m_{pK#pi}^{DCAFitter} (GeV/#it{c}^{2})", Variable::FromString("Candidates.fLiteM"), {240, 2.18, 2.42}}, {"m_{pK#pi}^{KF} (GeV/#it{c}^{2})", Variable::FromString("Candidates.fKFMassInv"), {240, 2.18, 2.42}});
}

int main(int argc, char* argv[]){
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./kf_vs_dca_qa filelistname (isMc=false)" << std::endl;
    exit(EXIT_FAILURE);
  }

  const bool isMc = argc > 2 ? HelperFunctions::StringToBool(argv[2]) : false;

  const std::string filelistname = argv[1];
  kf_vs_dca_qa(filelistname, isMc);

  return 0;
}