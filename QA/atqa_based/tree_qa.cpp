#include <string>

#include "AnalysisTree/TaskManager.hpp"
#include "AnalysisTree/Variable.hpp"
#include "Task.hpp"

using namespace AnalysisTree;

void tree_qa(const std::string& filelist){
  auto* man = TaskManager::GetInstance();

  auto* task = new QA::Task;
  task->SetOutputFileName("tree_qa.root");

  struct DataType {
    std::string name_;
    SimpleCut cut_;
  };

  std::vector<DataType> datatypes {
    {"background", EqualsCut("Candidates.KF_fSigBgStatus", 0)},
    {"prompt",     EqualsCut("Candidates.KF_fSigBgStatus", 1)},
    {"nonprompt",  EqualsCut("Candidates.KF_fSigBgStatus", 2)},
    {"wrongswap",  EqualsCut("Candidates.KF_fSigBgStatus", 3)},
    {"data",       EqualsCut("Candidates.KF_fSigBgStatus", -999)},
    {"impossible", SimpleCut({"Candidates.KF_fSigBgStatus"}, [](std::vector<double> par){ return par[0] != 0 && par[0] != 1 && par[0] != 2 && par[0] != 3 && par[0] != -999; })},
  };



  man->AddTask(task);
  man->Init({filelist}, {"aTree"});
  man->SetVerbosityPeriod(100);
  man->Run(-1);
  man->Finish();
}

int main(int argc, char* argv[]){
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./tree_qa filelistname" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string filelistname = argv[1];
  tree_qa(filelistname);

  return 0;
}
