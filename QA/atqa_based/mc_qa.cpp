#include <string>

#include "AnalysisTree/TaskManager.hpp"
#include "AnalysisTree/Variable.hpp"
#include "Task.hpp"

using namespace AnalysisTree;

template<typename T>
std::string to_string_with_precision(const T a_value, const int n);

void mc_qa(const std::string& filelist){
  auto* man = TaskManager::GetInstance();

  auto* task = new QA::Task;
  task->SetOutputFileName("mc_qa.root");

  const std::vector<SimpleCut> promptnesses {
    EqualsCut("Candidates.KF_fSigBgStatus", 1, "prompt"),
//     EqualsCut("Candidates.KF_fSigBgStatus", 2, "nonprompt")
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

  const std::vector<float> pTRanges{0.f, 2.f, 4.f, 6.f, 8.f, 12.f, 16.f};
  std::vector<SimpleCut> pTCuts;
  for(int iRange=0; iRange<pTRanges.size()-1; iRange++) {
    const std::string cutName = "pTsim_" + to_string_with_precision(pTRanges.at(iRange), 2) + "_" + to_string_with_precision(pTRanges.at(iRange+1), 2);
    pTCuts.emplace_back(RangeCut("Simulated.Sim_fPt", pTRanges.at(iRange), pTRanges.at(iRange+1), cutName));
  }

  const std::vector<float> pRanges{0.f, 2.f, 4.f, 6.f, 8.f, 12.f, 16.f};
  std::vector<SimpleCut> pCuts;
  for(int iRange=0; iRange<pRanges.size()-1; iRange++) {
    const std::string cutName = "psim_" + to_string_with_precision(pRanges.at(iRange), 2) + "_" + to_string_with_precision(pRanges.at(iRange+1), 2);
    pCuts.emplace_back(RangeCut("Simulated.Sim_fP", pRanges.at(iRange), pRanges.at(iRange+1), cutName));
  }

  const std::vector<float> lRanges{0.f, 0.02f, 0.04f, 0.06f, 0.08f, 0.12f, 0.16f, 0.20f};
  std::vector<SimpleCut> lCuts;
  for(int iRange=0; iRange<lRanges.size()-1; iRange++) {
    const std::string cutName = "lsim_" + to_string_with_precision(lRanges.at(iRange), 2) + "_" + to_string_with_precision(lRanges.at(iRange+1), 2);
    lCuts.emplace_back(RangeCut("Simulated.Sim_fDecayL", lRanges.at(iRange), lRanges.at(iRange+1), cutName));
  }

  const std::vector<float> tRanges{0.f, 0.2f, 0.4f, 0.6f, 0.8f, 1.2f, 1.6f, 2.f};
  std::vector<SimpleCut> tCuts;
  for(int iRange=0; iRange<tRanges.size()-1; iRange++) {
    const std::string cutName = "tsim_" + to_string_with_precision(tRanges.at(iRange), 2) + "_" + to_string_with_precision(tRanges.at(iRange+1), 2);
    tCuts.emplace_back(RangeCut("Simulated.Sim_fDecayT", tRanges.at(iRange), tRanges.at(iRange+1), cutName));
  }

  const int nbins = 400;
  std::vector<Quantity> vars {
    {"P",    "KF_fP",       "Sim_fP",       "KF_fDeltaP",  "p",      "GeV/c", {nbins, 0,    16 }, {nbins, -1, 1      }, {nbins, -5,  5 }, pCuts },
    {"Pt",   "KF_fPt",      "Sim_fPt",      "KF_fDeltaPt", "p_{T}",  "GeV/c", {nbins, 0,    16 }, {nbins, -1, 1      }, {nbins, -5,  5 }, pTCuts},
    {"Xsv",  "KF_fX",       "Sim_fDecayX",  "KF_fErrX",    "X_{SV}", "cm",    {nbins, -0.3, 0.3}, {nbins, -0.05, 0.05}, {nbins, -5,  5 }, pCuts },
    {"Ysv",  "KF_fY",       "Sim_fDecayY",  "KF_fErrY",    "Y_{SV}", "cm",    {nbins, -0.3, 0.3}, {nbins, -0.05, 0.05}, {nbins, -5,  5 }, pCuts },
    {"Zsv",  "KF_fZ",       "Sim_fDecayZ",  "KF_fErrZ",    "Z_{SV}", "cm",    {nbins, -20,  20 }, {nbins, -0.05, 0.05}, {nbins, -5,  5 }, pCuts },
    {"Xpv",  "Lite_fPosX",  "Sim_fEventX",  "KF_fErrPVX",  "X_{PV}", "cm",    {nbins, -0.3, 0.3}, {nbins, -0.05, 0.05}, {nbins, -5,  5 }, {}    },
    {"Ypv",  "Lite_fPosY",  "Sim_fEventY",  "KF_fErrPVY",  "Y_{PV}", "cm",    {nbins, -0.3, 0.3}, {nbins, -0.05, 0.05}, {nbins, -5,  5 }, {}    },
    {"Zpv",  "Lite_fPosZ",  "Sim_fEventZ",  "KF_fErrPVZ",  "Z_{PV}", "cm",    {nbins, -20,  20 }, {nbins, -0.05, 0.05}, {nbins, -5,  5 }, {}    },
    {"L",    "KF_fL",       "Sim_fDecayL",  "KF_fDeltaL",  "L",      "cm",    {nbins, -0.1, 0.2}, {nbins, -0.05, 0.05}, {nbins, -5,  5 }, lCuts },
    {"T",    "KF_fT",       "Sim_fDecayT",  "KF_fDeltaT",  "T",      "ps",    {nbins, -1,   2  }, {nbins, -2, 2      }, {nbins, -10, 10}, tCuts },
  };

  struct Histogram {
    std::string name_;
    std::string xaxis_title_;
    Variable variable_;
    TAxis xaxis_;
  };

  for(auto& pro : promptnesses) {
    for(auto& var : vars) {
      var.slice_cuts_.emplace_back(SimpleCut({"Simulated.ones"}, [](std::vector<double> par){return true;}, "total"));
      for(auto& slc : var.slice_cuts_) {
        Cuts* cutSelected = new Cuts((pro.GetTitle() + "_" + slc.GetTitle()).c_str(), {simpleCutSelected, pro, slc});
        Variable varResidual(("res_" + var.name_).c_str(), {{"Candidates", var.name_in_tree_rec_.c_str()}, {"Simulated", var.name_in_tree_mc_}}, []( std::vector<double>& v ) { return v.at(0) - v.at(1); });
        Variable varPull(("pull_" + var.name_).c_str(), {{"Candidates", var.name_in_tree_rec_.c_str()}, {"Simulated", var.name_in_tree_mc_}, {"Candidates", var.name_in_tree_error_.c_str()}}, []( std::vector<double>& v ) { return (v.at(0) - v.at(1)) / v.at(2); });
        std::vector<Histogram> histos {
          {("mc_" + var.name_ + "_" + slc.GetTitle()).c_str(), (var.title_ + "^{mc}, " + var.unit_).c_str(), Variable::FromString(("Simulated." + var.name_in_tree_mc_).c_str()), var.axis_plain_},
          {("rec_" + var.name_ + "_" + slc.GetTitle()).c_str(), (var.title_ + "^{rec}, " + var.unit_).c_str(), Variable::FromString(("Candidates." + var.name_in_tree_rec_).c_str()), var.axis_plain_},
          {("res_" + var.name_ + "_" + slc.GetTitle()).c_str(), (var.title_ + "^{rec} - " + var.title_ + "^{mc}, " + var.unit_).c_str(), varResidual, var.axis_res_},
          {("pull_" + var.name_ + "_" + slc.GetTitle()).c_str(), ("(" + var.title_ + "^{rec} - " + var.title_ + "^{mc}) / #sigma_{" + var.title_ + "^{rec}}").c_str(), varPull, var.axis_pull_},
        };

        for(auto& histo : histos) {
          task->AddH1(histo.name_, {histo.xaxis_title_, histo.variable_, histo.xaxis_}, cutSelected);
        }

        if(histos.at(0).name_ == ("mc_" + var.name_).c_str() && histos.at(1).name_ == ("rec_" + var.name_).c_str()) {
          task->AddH2(("corr_" + var.name_).c_str(), {histos.at(0).xaxis_title_, histos.at(0).variable_, var.axis_plain_}, {histos.at(1).xaxis_title_, histos.at(1).variable_, var.axis_plain_}, cutSelected);
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

template<typename T>
std::string to_string_with_precision(const T a_value, const int n) {
  std::ostringstream out;
  out.precision(n);
  out << std::fixed << a_value;
  return out.str();
}
