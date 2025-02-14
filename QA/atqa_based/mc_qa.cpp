#include <string>

#include "AnalysisTree/Constants.hpp"
#include "AnalysisTree/HelperFunctions.hpp"
#include "AnalysisTree/TaskManager.hpp"
#include "AnalysisTree/Variable.hpp"
#include "Task.hpp"

using namespace AnalysisTree;

void mc_qa(const std::string& filelist){
  auto* man = TaskManager::GetInstance();

  auto* task = new QA::Task;
  task->SetOutputFileName("mc_qa.root");

  const std::vector<SimpleCut> promptnesses {
    EqualsCut("Candidates.KF_fSigBgStatus", 1, "prompt"),
    EqualsCut("Candidates.KF_fSigBgStatus", 2, "nonprompt")
  };

//   SimpleCut simpleCutSelected = RangeCut("Candidates.KF_fIsSelected", 0.9, 1.1);
  const SimpleCut simpleCutSelected = RangeCut("Candidates.KF_fIsSelected", -0.1, 1.1);

  struct Quantity {
    std::string name_;
    std::string name_in_tree_rec_;
    std::string name_in_tree_mc_;
    std::string name_in_tree_error_;
    std::string title_;
    std::string unit_;
    TAxis axis_plain_;
    TAxis axis_res_;
    TAxis axis_pull_;
    std::vector<SimpleCut> slice_cuts_;
  };

  auto pTCuts = HelperFunctions::CreateRangeCuts({0.f, 2.f, 4.f, 6.f, 8.f, 12.f, 16.f}, "pTsim_", "Simulated.Sim_fPt", 2);
  auto pCuts = HelperFunctions::CreateRangeCuts({0.f, 2.f, 4.f, 6.f, 8.f, 12.f, 16.f}, "psim_", "Simulated.Sim_fP", 2);
  auto lCuts = HelperFunctions::CreateRangeCuts({0.f, 0.02f, 0.04f, 0.06f, 0.08f, 0.12f, 0.16f, 0.20f}, "lsim_", "Simulated.Sim_fLDecay", 2);
  auto tCuts = HelperFunctions::CreateRangeCuts({0.f, 0.2f, 0.4f, 0.6f, 0.8f, 1.2f, 1.6f, 2.f}, "tsim_", "Simulated.Sim_fTDecay", 2);
  auto nPVCCuts = HelperFunctions::CreateEqualCuts({0, 1, 2, 3}, "nPVC_", "Candidates.Lite_fNProngsContributorsPV", 0);

  const int nbins = 400;

  std::vector<Quantity> vars {
    {"P",    "KF_fP",           "Sim_fP",      "KF_fErrP",             "p",      "GeV/c", {nbins, 0,    16 }, {nbins, -1, 1      }, {nbins, -5,  5 }, pCuts   },
    {"Pt",   "KF_fPt",          "Sim_fPt",     "KF_fErrPt",            "p_{T}",  "GeV/c", {nbins, 0,    16 }, {nbins, -1, 1      }, {nbins, -5,  5 }, pTCuts  },
    {"Xsv",  "KF_fX",           "Sim_fXDecay", "KF_fErrX",             "X_{SV}", "cm",    {nbins, -0.5, 0.5}, {nbins, -0.05, 0.05}, {nbins, -5,  5 }, pCuts   },
    {"Ysv",  "KF_fY",           "Sim_fYDecay", "KF_fErrY",             "Y_{SV}", "cm",    {nbins, -0.5, 0.5}, {nbins, -0.05, 0.05}, {nbins, -5,  5 }, pCuts   },
    {"Zsv",  "KF_fZ",           "Sim_fZDecay", "KF_fErrZ",             "Z_{SV}", "cm",    {nbins, -20,  20 }, {nbins, -0.05, 0.05}, {nbins, -5,  5 }, pCuts   },
    {"Xpv",  "Lite_fPosX",      "Sim_fXEvent", "KF_fErrPVX",           "X_{PV}", "cm",    {nbins, -0.1, 0.1}, {nbins, -0.05, 0.05}, {nbins, -5,  5 }, nPVCCuts},
    {"Ypv",  "Lite_fPosY",      "Sim_fYEvent", "KF_fErrPVY",           "Y_{PV}", "cm",    {nbins, -0.1, 0.1}, {nbins, -0.05, 0.05}, {nbins, -5,  5 }, nPVCCuts},
    {"Zpv",  "Lite_fPosZ",      "Sim_fZEvent", "KF_fErrPVZ",           "Z_{PV}", "cm",    {nbins, -20,  20 }, {nbins, -0.05, 0.05}, {nbins, -5,  5 }, nPVCCuts},
    {"L",    "KF_fDecayLength", "Sim_fLDecay", "KF_fDecayLengthError", "L",      "cm",    {nbins, -0.1, 0.2}, {nbins, -0.05, 0.05}, {nbins, -5,  5 }, lCuts   },
    {"T",    "KF_fT",           "Sim_fTDecay", "KF_fErrT",             "T",      "ps",    {nbins, -1,   2  }, {nbins, -2, 2      }, {nbins, -10, 10}, tCuts   },
  };


  struct Histogram {
    std::string name_;
    std::string xaxis_title_;
    Variable variable_;
    TAxis xaxis_;
    bool to_be_sliced_;
  };

  for(auto& pro : promptnesses) {
    Cuts* cutsTotal = new Cuts(pro.GetTitle() + "_total", {simpleCutSelected, pro});
    for(auto& var : vars) {

      Variable varMc(Variable::FromString("Simulated." + var.name_in_tree_mc_));
      Variable varRec(Variable::FromString("Candidates." + var.name_in_tree_rec_));
      Variable varResidual("res_" + var.name_, {{"Candidates", var.name_in_tree_rec_}, {"Simulated", var.name_in_tree_mc_}}, []( std::vector<double>& v ) { return v.at(0) - v.at(1); });
      Variable varPull("pull_" + var.name_, {{"Candidates", var.name_in_tree_rec_}, {"Simulated", var.name_in_tree_mc_}, {"Candidates", var.name_in_tree_error_}}, []( std::vector<double>& v ) { return (v.at(0) - v.at(1)) / v.at(2); });

      std::vector<Histogram> histos {
        {"mc_" + var.name_, var.title_ + "^{mc}, " + var.unit_, varMc, var.axis_plain_, false},
        {"rec_" + var.name_, var.title_ + "^{rec}, " + var.unit_, varRec, var.axis_plain_, false},
        {"res_" + var.name_, var.title_ + "^{rec} - " + var.title_ + "^{mc}, " + var.unit_, varResidual, var.axis_res_, true},
        {"pull_" + var.name_, "(" + var.title_ + "^{rec} - " + var.title_ + "^{mc}) / #sigma_{" + var.title_ + "^{rec}}", varPull, var.axis_pull_, true},
      };

      for(auto& histo : histos) {
        task->AddH1(histo.name_, {histo.xaxis_title_, histo.variable_, histo.xaxis_}, cutsTotal);
      }

      if(histos.at(0).name_ == "mc_" + var.name_ && histos.at(1).name_ == "rec_" + var.name_) {
        task->AddH2("corr_" + var.name_, {histos.at(0).xaxis_title_, histos.at(0).variable_, var.axis_plain_}, {histos.at(1).xaxis_title_, histos.at(1).variable_, var.axis_plain_}, cutsTotal);
      }

      for(auto& slc : var.slice_cuts_) {
        Cuts* cutSlice = new Cuts(pro.GetTitle() + "_" + slc.GetTitle(), {simpleCutSelected, pro, slc});

        for(auto& histo : histos) {
          if(histo.name_ == "mc_" + var.name_ || histo.name_ == "rec_" + var.name_) continue;
          task->AddH1(histo.name_ + "_" + slc.GetTitle(), {histo.xaxis_title_, histo.variable_, histo.xaxis_}, cutSlice);
        }
      } // var.slice_cuts_
    } // vars
  } // promptnesses

  man->AddTask(task);
  man->Init({filelist}, {"aTree"});
  man->SetVerbosityPeriod(100);
  man->Run(-1);
  man->Finish();
}

int main(int argc, char* argv[]){
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./mc_qa filelistname" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string filelistname = argv[1];
  mc_qa(filelistname);

  return 0;
}

