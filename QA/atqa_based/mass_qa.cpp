//
// Created by oleksii on 18.03.25.
//

#include "Task.hpp"

#include "AnalysisTree/HelperFunctions.hpp"
#include "AnalysisTree/TaskManager.hpp"

using namespace AnalysisTree;

void MassQA(QA::Task& task, TH1D* hEff);

double Interpolate(const TH1* h, double value);

std::vector<SimpleCut> datatypes {
    RangeCut("Candidates.KF_fSigBgStatus", -0.1, 3.1, "all"),
    EqualsCut("Candidates.KF_fSigBgStatus", 0,"background"),
    EqualsCut("Candidates.KF_fSigBgStatus", 1,"prompt"),
    EqualsCut("Candidates.KF_fSigBgStatus", 2,"nonprompt"),
    EqualsCut("Candidates.KF_fSigBgStatus", 3,"wrongswap"),
    //   EqualsCut("Candidates.KF_fSigBgStatus", -999, "data"),
    //   SimpleCut({"Candidates.KF_fSigBgStatus"}, [](std::vector<double> par){ return par[0] != 0 && par[0] != 1 && par[0] != 2 && par[0] != 3 && par[0] != -999; }, "impossible"),
};

void mass_qa(const std::string& filelistname, const std::string& fileEffMapName, const std::string& histoEffMapName) {
  auto* man = TaskManager::GetInstance();

  TFile* fileEff = TFile::Open(fileEffMapName.c_str(), "open");
  TH1D* hEff = fileEff->Get<TH1D>(histoEffMapName.c_str());

  auto* task = new QA::Task;
  task->SetOutputFileName("mass_qa.root");

  MassQA(*task, hEff);

  man->AddTask(task);
  man->Init({filelistname}, {"aTree"});
  man->SetVerbosityPeriod(1000);
  man->Run();
  man->Finish();

  fileEff->Close();
}

void MassQA(QA::Task& task, TH1D* hEff) {
  auto TCuts = HelperFunctions::CreateRangeCuts({0.2, 0.35, 0.5, 0.7, 0.9, 1.6}, "T_", "Candidates.KF_fT", 2);

  const std::string varName = "Mass";
  const std::string axisTitle = "m_{pK#pi} (GeV/#it{c}^{2})";
  const std::string varNameInTree = "Candidates.KF_fMassInv";
  const TAxis varAxis = {600, 1.98, 2.58};

  Variable recEffWeigth("recEffWeigth",
                        {{"Candidates", "KF_fT"}},
                        [=] (std::vector<double>& par) {
                          const double eff = Interpolate(hEff, par[0]);
                          return eff < 1e-4 ? 0. : 1. / eff;
                        } );

  for(auto& dt : datatypes) {
    task.SetTopLevelDirName(dt.GetTitle());
    Cuts* cutsTotal = new Cuts(dt.GetTitle() + "_total", {dt});
    task.AddH1(varName + "_" + dt.GetTitle(), {axisTitle, Variable::FromString(varNameInTree), varAxis}, cutsTotal);
    task.AddH1(varName + "_W_" + dt.GetTitle(), {axisTitle, Variable::FromString(varNameInTree), varAxis}, cutsTotal, recEffWeigth);

    for(auto& slc : TCuts) {
      Cuts* cutSlice = new Cuts(dt.GetTitle() + "_" + slc.GetTitle(), {dt, slc});
      task.AddH1(varName + "_" + dt.GetTitle() + "_" + slc.GetTitle(), {axisTitle, Variable::FromString(varNameInTree), varAxis}, cutSlice);
      task.AddH1(varName + "_W_" + dt.GetTitle() + "_" + slc.GetTitle(), {axisTitle, Variable::FromString(varNameInTree), varAxis}, cutSlice, recEffWeigth);
    } // TCuts
  } // datatypes

}

double Interpolate(const TH1* h, double value) {
  double result;
  if(value <= h->GetBinLowEdge(1) || value >= h->GetBinLowEdge(h->GetNbinsX()+1)) result = 0.;
  else result = h->Interpolate(value);

  return result;
}

int main(int argc, char* argv[]){
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./mass_qa filelistname, fileEffMapName, histoEffMapName" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string filelistname = argv[1];
  const std::string fileEffMapName = argv[2];
  const std::string histoEffMapName = argv[3];

  mass_qa(filelistname, fileEffMapName, histoEffMapName);

  return 0;
}