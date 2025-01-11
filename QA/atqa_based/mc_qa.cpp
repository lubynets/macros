#include <string>

#include "AnalysisTree/TaskManager.hpp"
#include "AnalysisTree/Variable.hpp"
#include "Task.hpp"

using namespace AnalysisTree;

void example(const std::string& filelist){
  auto* man = TaskManager::GetInstance();

  auto* task = new QA::Task;
  task->SetOutputFileName("mc_qa.root");

  struct Promptness {
    std::string name_;
    SimpleCut cut_;
  };

  std::vector<Promptness> promptnesses {
    {"prompt",    EqualsCut("Candidates.KF_fSigBgStatus", 1)},
    {"nonprompt", EqualsCut("Candidates.KF_fSigBgStatus", 2)}
  };

  SimpleCut simpleCutSelected = EqualsCut("Candidates.KF_fIsSelected", 1);

  struct Quantity {
    std::string name_;
    std::string name_in_tree_rec_;
    std::string name_in_tree_mc_;
    std::string name_in_tree_error_;
    std::string title_;
    std::string unit_;
    int nbins_plain_;
    float xlow_plain_;
    float xup_plain_;
    int nbins_res_;
    float xlow_res_;
    float xup_res_;
    int nbins_pull_;
    float xlow_pull_;
    float xup_pull_;
  };

  const int nbins = 400;

  std::vector<Quantity> vars {
    {"P",   "KF_fP",      "Sim_fP",       "KF_fDeltaP",  "p",      "GeV/c", nbins, 0,    12,  nbins, -1, 1,       nbins, -5,  5 },
    {"Pt",  "KF_fPt",     "Sim_fPt",      "KF_fDeltaPt", "p_{T}",  "GeV/c", nbins, 0,    12,  nbins, -1, 1,       nbins, -5,  5 },
    {"Xsv", "KF_fX",      "Sim_fDecayX",  "KF_fErrX",    "X_{SV}", "cm",    nbins, -0.3, 0.3, nbins, -0.05, 0.05, nbins, -5,  5 },
    {"Ysv", "KF_fY",      "Sim_fDecayY",  "KF_fErrY",    "Y_{SV}", "cm",    nbins, -0.3, 0.3, nbins, -0.05, 0.05, nbins, -5,  5 },
    {"Zsv", "KF_fZ",      "Sim_fDecayZ",  "KF_fErrZ",    "Z_{SV}", "cm",    nbins, -20,  20,  nbins, -0.05, 0.05, nbins, -5,  5 },
    {"Xpv", "Lite_fPosX", "Sim_fEventX",  "KF_fErrPVX",  "X_{PV}", "cm",    nbins, -0.3, 0.3, nbins, -0.05, 0.05, nbins, -5,  5 },
    {"Ypv", "Lite_fPosY", "Sim_fEventY",  "KF_fErrPVY",  "Y_{PV}", "cm",    nbins, -0.3, 0.3, nbins, -0.05, 0.05, nbins, -5,  5 },
    {"Zpv", "Lite_fPosZ", "Sim_fEventZ",  "KF_fErrPVZ",  "Z_{PV}", "cm",    nbins, -20,  20,  nbins, -0.05, 0.05, nbins, -5,  5 },
    {"L",   "KF_fL",      "Sim_fDecayL",  "KF_fDeltaL",  "L",      "cm",    nbins, -0.1, 0.2, nbins, -0.05, 0.05, nbins, -5,  5 },
    {"T",   "KF_fT",      "Sim_fDecayT",  "KF_fDeltaT",  "T",      "ps",    nbins, -1,   2,   nbins, -2, 2,       nbins, -10, 10},
  };

  for(auto& pro : promptnesses) {
    Cuts* cutSelected = new Cuts(pro.name_, {/*simpleCutSelected,*/ pro.cut_});
    for(auto& var : vars) {
      task->AddH1({(var.title_ + "^{mc}, " + var.unit_).c_str(), Variable::FromString(("Simulated." + var.name_in_tree_mc_).c_str()), {nbins, var.xlow_plain_, var.xup_plain_}}, cutSelected);
      task->AddH1({(var.title_ + "^{rec}, " + var.unit_).c_str(), Variable::FromString(("Candidates." + var.name_in_tree_rec_).c_str()), {nbins, var.xlow_plain_, var.xup_plain_}}, cutSelected);
      task->AddH2({(var.title_ + "^{mc}, " + var.unit_).c_str(), Variable::FromString(("Simulated." + var.name_in_tree_mc_).c_str()), {nbins, var.xlow_plain_, var.xup_plain_}}, {(var.title_ + "^{rec}, " + var.unit_).c_str(), Variable::FromString(("Candidates." + var.name_in_tree_rec_).c_str()), {nbins, var.xlow_plain_, var.xup_plain_}}, cutSelected);

      Variable residual(("res_" + var.name_).c_str(), {{"Candidates", var.name_in_tree_rec_.c_str()}, {"Simulated", var.name_in_tree_mc_}}, []( std::vector<double>& v ) { return v.at(0) - v.at(1); });
      Variable pull(("pull_" + var.name_).c_str(), {{"Candidates", var.name_in_tree_rec_.c_str()}, {"Simulated", var.name_in_tree_mc_}, {"Candidates", var.name_in_tree_error_.c_str()}}, []( std::vector<double>& v ) { return (v.at(0) - v.at(1)) / v.at(2); });

      task->AddH1({(var.title_ + "^{rec} - " + var.title_ + "^{mc}, " + var.unit_).c_str(), residual, {nbins, var.xlow_res_, var.xup_res_}}, cutSelected);
      task->AddH1({("(" + var.title_ + "^{rec} - " + var.title_ + "^{mc}) / #sigma_{" + var.title_ + "^{rec}}").c_str(), pull, {nbins, var.xlow_pull_, var.xup_pull_}}, cutSelected);
    }
  }

  man->AddTask(task);
  man->Init({filelist}, {"aTree"});
  man->SetVerbosityPeriod(100);
  man->Run(-1);
  man->Finish();
}

int main(int argc, char* argv[]){
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./example filelistname" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string filelistname = argv[1];
  example(filelistname);

  return 0;
}
