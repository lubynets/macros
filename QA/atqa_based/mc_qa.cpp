#include <string>

#include "AnalysisTree/Constants.hpp"
#include "AnalysisTree/HelperFunctions.hpp"
#include "AnalysisTree/TaskManager.hpp"
#include "AnalysisTree/Variable.hpp"
#include "Task.hpp"

using namespace AnalysisTree;

void PullsAndResidualsQA(QA::Task& task);
void EfficiencyQA(QA::Task& task);

std::vector selectCuts {
  RangeCut("Candidates.KF_fIsSelected", 0.9, 1.1, "isSel"),
  RangeCut("Candidates.KF_fIsSelected", -0.1, 1.1, "noSel")
};

void mc_qa(const std::string& filelist){
  auto* man = TaskManager::GetInstance();

  auto* task = new QA::Task;
  task->SetOutputFileName("mc_qa.root");

  PullsAndResidualsQA(*task);
  EfficiencyQA(*task);

  man->AddTask(task);
  man->Init({filelist}, {"aTree"});
  man->SetVerbosityPeriod(10000);
  man->Run();
  man->Finish();
}

void EfficiencyQA(QA::Task& task) {
  task.SetTopLevelDirName("EfficiencyQA");

  struct Quantity {
    std::string name_;
    std::string name_in_tree_rec_;
    std::string name_in_tree_gen_;
    std::string title_;
    std::string unit_;
    TAxis axis_;
    std::vector<SimpleCut> cut_rec_;
    std::vector<SimpleCut> cut_gen_;
  };

  const int nbins = 400;
  const std::pair<float, float> rapidityCut = std::make_pair(-0.8, 0.8);
  SimpleCut rapidityCutRec = RangeCut(Variable::FromString("Candidates.Lite_fY"), rapidityCut.first, rapidityCut.second);
  SimpleCut rapidityCutGen = RangeCut(Variable::FromString("Generated.Gen_fY"), rapidityCut.first, rapidityCut.second);
  const std::vector<Quantity> vars {
    {"Y",  "Lite_fY",         "Gen_fY",      "y",     "",      {nbins, -1, 1}, {},               {}              },
    {"Pt", "KF_fPt",          "Gen_fPt",     "p_{T}", "GeV/c", {nbins, 0, 16}, {rapidityCutRec}, {rapidityCutGen}},
    {"L",  "KF_fDecayLength", "Gen_fLDecay", "L",     "cm",    {nbins, 0, 2},  {rapidityCutRec}, {rapidityCutGen}},
    {"T",  "KF_fT",           "Gen_fTDecay", "T",     "ps",    {nbins, 0, 2},  {rapidityCutRec}, {rapidityCutGen}},
  };

  struct Promptness {
    std::string name_;
    SimpleCut cut_rec_;
    SimpleCut cut_gen_;
  };

  const std::vector<Promptness> promptnesses{
    {"prompt", EqualsCut("Candidates.KF_fSigBgStatus", 1), EqualsCut("Generated.Gen_fOriginMcGen", 1)},
    {"nonprompt", EqualsCut("Candidates.KF_fSigBgStatus", 2), EqualsCut("Generated.Gen_fOriginMcGen", 2)},
  };

  for (auto& promptness : promptnesses) {
    task.SetTopLevelDirName("EfficiencyQA/" + promptness.name_ + "/gen");
    Cuts* genCut2D = new Cuts(promptness.name_, {promptness.cut_gen_});
    task.AddH2("gen_" + vars.at(0).name_ + "_Vs_" + vars.at(1).name_,
              {vars.at(0).title_ + vars.at(0).unit_, Variable::FromString("Generated." + vars.at(0).name_in_tree_gen_), vars.at(0).axis_},
              {vars.at(1).title_ + ", " + vars.at(1).unit_, Variable::FromString("Generated." + vars.at(1).name_in_tree_gen_), vars.at(1).axis_}, genCut2D);
    for (auto& var : vars) {
      task.SetTopLevelDirName("EfficiencyQA/" + promptness.name_ + "/gen");
      const std::string comma = var.unit_ == "" ? "" : ", ";
      Cuts* genCut1D = new Cuts(promptness.name_, {promptness.cut_gen_});
      genCut1D->AddCuts(var.cut_gen_);
      task.AddH1("gen_" + var.name_, {var.title_ + comma + var.unit_, Variable::FromString("Generated." + var.name_in_tree_gen_), var.axis_}, genCut1D);
      for(auto& sel : selectCuts) {
        task.SetTopLevelDirName("EfficiencyQA/" + promptness.name_ + "/rec/" + sel.GetTitle());
        Cuts* recCut1D = new Cuts(promptness.name_, {promptness.cut_rec_, sel});
        recCut1D->AddCuts(var.cut_rec_);
        task.AddH1("rec_" + var.name_, {var.title_ + comma + var.unit_, Variable::FromString("Candidates." + var.name_in_tree_rec_), var.axis_}, recCut1D);
      } // selectCuts
    } // vars
    for(auto& sel : selectCuts) {
      task.SetTopLevelDirName("EfficiencyQA/" + promptness.name_ + "/rec/" + sel.GetTitle());
      Cuts* recCut2D = new Cuts(promptness.name_, {promptness.cut_rec_, sel});
      task.AddH2("rec_" + vars.at(0).name_ + "_Vs_" + vars.at(1).name_,
                {vars.at(0).title_ + vars.at(0).unit_, Variable::FromString("Candidates." + vars.at(0).name_in_tree_rec_), vars.at(0).axis_},
                {vars.at(1).title_ + ", " + vars.at(1).unit_, Variable::FromString("Candidates." + vars.at(1).name_in_tree_rec_), vars.at(1).axis_}, recCut2D);
    } // selectCuts
  } // promptnesses
} // EfficiencyQA()

void PullsAndResidualsQA(QA::Task& task) {
    const std::vector<SimpleCut> promptnesses {
    EqualsCut("Candidates.KF_fSigBgStatus", 1, "prompt"),
    EqualsCut("Candidates.KF_fSigBgStatus", 2, "nonprompt")
  };

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

  const int nbins = 400;
  std::vector<Quantity> varsCand {
    {"P",    "KF_fP",           "Sim_fP",      "KF_fErrP",             "p",      "GeV/c", {nbins, 0,    16 }, {nbins, -1, 1      }, {nbins, -5,  5 }, pCuts   },
    {"Pt",   "KF_fPt",          "Sim_fPt",     "KF_fErrPt",            "p_{T}",  "GeV/c", {nbins, 0,    16 }, {nbins, -1, 1      }, {nbins, -5,  5 }, pTCuts  },
    {"Xsv",  "KF_fX",           "Sim_fXDecay", "KF_fErrX",             "X_{SV}", "cm",    {nbins, -0.5, 0.5}, {nbins, -0.05, 0.05}, {nbins, -5,  5 }, pTCuts  },
    {"Ysv",  "KF_fY",           "Sim_fYDecay", "KF_fErrY",             "Y_{SV}", "cm",    {nbins, -0.5, 0.5}, {nbins, -0.05, 0.05}, {nbins, -5,  5 }, pTCuts  },
    {"Zsv",  "KF_fZ",           "Sim_fZDecay", "KF_fErrZ",             "Z_{SV}", "cm",    {nbins, -20,  20 }, {nbins, -0.05, 0.05}, {nbins, -5,  5 }, pTCuts  },
    {"L",    "KF_fDecayLength", "Sim_fLDecay", "KF_fDecayLengthError", "L",      "cm",    {nbins, -0.1, 0.2}, {nbins, -0.05, 0.05}, {nbins, -5,  5 }, lCuts   },
    {"T",    "KF_fT",           "Sim_fTDecay", "KF_fErrT",             "T",      "ps",    {nbins, -1,   2  }, {nbins, -2, 2      }, {nbins, -10, 10}, tCuts   },
  };

  constexpr float Shift001 = 0.01;
  auto nPVCCuts = HelperFunctions::CreateRangeCuts({0, 10-Shift001, 20-Shift001, 40-Shift001, 60-Shift001, 80-Shift001, 100}, "nPVC_", "Events.Ev_fMultNTracksPV", 2);
  std::vector<Quantity> varsEve {
    {"Xpv", "Ev_fPosX", "Ev_fMcPosX", "Ev_fPosXErr", "X_{PV}", "cm", {nbins, -0.1, 0.1}, {nbins, -0.05, 0.05}, {nbins, -5,  5 }, nPVCCuts},
    {"Ypv", "Ev_fPosY", "Ev_fMcPosY", "Ev_fPosYErr", "Y_{PV}", "cm", {nbins, -0.1, 0.1}, {nbins, -0.05, 0.05}, {nbins, -5,  5 }, nPVCCuts},
    {"Zpv", "Ev_fPosZ", "Ev_fMcPosZ", "Ev_fPosZErr", "Z_{PV}", "cm", {nbins, -20,  20 }, {nbins, -0.05, 0.05}, {nbins, -5,  5 }, nPVCCuts},
  };

  struct Histogram {
    std::string name_;
    std::string xaxis_title_;
    Variable variable_;
    TAxis xaxis_;
    bool to_be_sliced_;
  };

  auto BuildPullsAndResiduals = [&] (Quantity& var,
                                    const std::string& recTreeName,
                                    const std::string& mcTreeName,
                                    const std::string& errTreeName,
                                    Cuts* cutExternal = nullptr) {
    Variable varMc(Variable::FromString(mcTreeName + "." + var.name_in_tree_mc_));
    Variable varRec(Variable::FromString(recTreeName + "." + var.name_in_tree_rec_));
    Variable varResidual("res_" + var.name_, {{recTreeName, var.name_in_tree_rec_}, {mcTreeName, var.name_in_tree_mc_}}, []( std::vector<double>& v ) { return v.at(0) - v.at(1); });
    Variable varPull("pull_" + var.name_, {{recTreeName, var.name_in_tree_rec_}, {mcTreeName, var.name_in_tree_mc_}, {errTreeName, var.name_in_tree_error_}}, []( std::vector<double>& v ) { return (v.at(0) - v.at(1)) / v.at(2); });

    std::vector<Histogram> histos {
      {"mc_" + var.name_, var.title_ + "^{mc}, " + var.unit_, varMc, var.axis_plain_, false},
      {"rec_" + var.name_, var.title_ + "^{rec}, " + var.unit_, varRec, var.axis_plain_, false},
      {"res_" + var.name_, var.title_ + "^{rec} - " + var.title_ + "^{mc}, " + var.unit_, varResidual, var.axis_res_, true},
      {"pull_" + var.name_, "(" + var.title_ + "^{rec} - " + var.title_ + "^{mc}) / #sigma_{" + var.title_ + "^{rec}}", varPull, var.axis_pull_, true},
    };

    for(auto& histo : histos) {
      task.AddH1(histo.name_, {histo.xaxis_title_, histo.variable_, histo.xaxis_}, cutExternal);
    }

    if(histos.at(0).name_ == "mc_" + var.name_ && histos.at(1).name_ == "rec_" + var.name_) {
      task.AddH2("corr_" + var.name_, {histos.at(0).xaxis_title_, histos.at(0).variable_, var.axis_plain_}, {histos.at(1).xaxis_title_, histos.at(1).variable_, var.axis_plain_}, cutExternal);
    }

    for(auto& slc : var.slice_cuts_) {
      const std::string cutExternalName = cutExternal != nullptr ? cutExternal->GetName() + "_" : "";
      Cuts* cutSlice = new Cuts(cutExternalName + slc.GetTitle(), {slc});
      if(cutExternal != nullptr) cutSlice->AddCuts(cutExternal->GetCuts());

      for(auto& histo : histos) {
        if(histo.name_ == "mc_" + var.name_ || histo.name_ == "rec_" + var.name_) continue;
        task.AddH1(histo.name_ + "_" + slc.GetTitle(), {histo.xaxis_title_, histo.variable_, histo.xaxis_}, cutSlice);
      }
    } // var.slice_cuts_
  };

  for(auto& sel : selectCuts) {
    for(auto& pro : promptnesses) {
      task.SetTopLevelDirName("PullsAndResiduals/CandidatesQA/" + sel.GetTitle() + "/" + pro.GetTitle());
      Cuts* cutsTotal = new Cuts(pro.GetTitle() + "_total", {sel, pro});
      for(auto& var : varsCand) {
        BuildPullsAndResiduals(var, "Candidates", "Simulated", "Candidates", cutsTotal);
      } // varsCand
    } // promptnesses
  } // selectCuts

  task.SetTopLevelDirName("PullsAndResiduals/EventsQA");
  for(auto& var : varsEve) {
    BuildPullsAndResiduals(var, "Events", "Events", "Events");
  } // varsEve
} // PullsAndResidualsQA()

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

