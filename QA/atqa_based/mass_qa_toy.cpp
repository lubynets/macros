//
// Created by oleksii on 29.08.2025.
//
#include "Task.hpp"

#include "AnalysisTree/HelperFunctions.hpp"
#include "AnalysisTree/TaskManager.hpp"

using namespace AnalysisTree;

const TAxis massAxis = {600, 1.98, 2.58};
const std::string massAxisTitle = "m_{pK#pi} (GeV/#it{c}^{2})";

const TAxis pTAxis = {200, 0, 20};
const std::string pTAxisTitle = "#it{p}_{T} (GeV/#it{c})";

void MassQABdtToy(QA::Task& task) {
  task.SetTopLevelDirName("");

  task.AddH1("hPt", {pTAxisTitle, Variable::FromString("Candidates.fLitePt"), pTAxis});

  for(int iPt=9; iPt<10; ++iPt) {
    const float lo = 0.1*iPt;
    const float hi = 0.1*(iPt+1);
    SimpleCut pTCut = RangeCut("Candidates.fLitePt", lo, hi);
    Cuts* cuts = new Cuts("cuts", {pTCut});
    const std::string histoName = "hMass_" + HelperFunctions::ToStringWithPrecision(lo, 1) + "_" + HelperFunctions::ToStringWithPrecision(hi, 1);
    task.AddH1(histoName, {massAxisTitle, Variable::FromString("Candidates.fKFMassInv"), massAxis}, cuts);
  }
}

void mass_qa_toy(const std::string& filelistname) {

  auto* man = TaskManager::GetInstance();
  const std::string fileOutName = "mass_qa_toy.root";

  auto* task = new QA::Task;
  task->SetOutputFileName(fileOutName);

  MassQABdtToy(*task);

  man->AddTask(task);
  man->Init({filelistname}, {"aTree"});
  man->SetVerbosityFrequency(10);
  man->Run();
  man->Finish();
}

int main(int argc, char* argv[]){
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./mass_qa_toy filelistname" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string filelistname = argv[1];

  mass_qa_toy(filelistname);

  return 0;
}