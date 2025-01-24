#include <string>

#include "AnalysisTree/Constants.hpp"
#include "AnalysisTree/HelperFunctions.hpp"
#include "AnalysisTree/TaskManager.hpp"
#include "AnalysisTree/Variable.hpp"
#include "Task.hpp"

using namespace AnalysisTree;

void PidQA(QA::Task& task);
void TopoQA(QA::Task& task);

std::vector<SimpleCut> datatypes {
//   EqualsCut("Candidates.KF_fSigBgStatus", 0,    "background"),
//   EqualsCut("Candidates.KF_fSigBgStatus", 1,    "prompt"),
//   EqualsCut("Candidates.KF_fSigBgStatus", 2,    "nonprompt"),
//   EqualsCut("Candidates.KF_fSigBgStatus", 3,    "wrongswap"),
  EqualsCut("Candidates.KF_fSigBgStatus", -999, "data"),
//   SimpleCut({"Candidates.KF_fSigBgStatus"}, [](std::vector<double> par){ return par[0] != 0 && par[0] != 1 && par[0] != 2 && par[0] != 3 && par[0] != -999; }, "impossible"),
};

std::vector<SimpleCut> dcaFitterSelections {
  RangeCut("Candidates.KF_fIsSelected", 0.9, 1.1, "isDcaFSel"),
  RangeCut("Candidates.KF_fIsSelected", -0.1, 1.1, "noDcaFSel"),
};

std::vector<SimpleCut> topoSelectionCuts{
  RangeCut("Candidates.KF_fChi2PrimProton",    5,          HugeNumber),
  RangeCut("Candidates.KF_fChi2PrimKaon",      5,          HugeNumber),
  RangeCut("Candidates.KF_fChi2PrimPion",      5,          HugeNumber),
  RangeCut("Candidates.KF_fChi2geoPionKaon",   -HugeNumber, 3        ),
  RangeCut("Candidates.KF_fChi2geoProtonKaon", -HugeNumber, 3        ),
  RangeCut("Candidates.KF_fChi2geoProtonPion", -HugeNumber, 3        ),
  RangeCut("Candidates.KF_fDCAPionKaon",       -HugeNumber, 0.01     ),
  RangeCut("Candidates.KF_fDCAProtonKaon",     -HugeNumber, 0.01     ),
  RangeCut("Candidates.KF_fDCAProtonPion",     -HugeNumber, 0.01     ),
  RangeCut("Candidates.KF_fChi2geo",           -HugeNumber, 3        ),
  RangeCut("Candidates.KF_fChi2topo",          -HugeNumber, 5        ),
  RangeCut("Candidates.KF_fLdL",               3,          HugeNumber),
};

void treeKF_qa(const std::string& filelist){
  auto* man = TaskManager::GetInstance();

  auto* task = new QA::Task;
  task->SetOutputFileName("treeKF_qa.root");

  TopoQA(*task);
  PidQA(*task);

  man->AddTask(task);
  man->Init({filelist}, {"aTree"});
  man->SetVerbosityPeriod(100);
  man->Run(-1);
  man->Finish();
}

void PidQA(QA::Task& task) {
  enum eProngSpecies : int {
    kProton = 0,
    kKaon,
    kPion,
    nProngSpecies
  };

  enum ePidDetectors : int {
    kTpc = 0,
    kTof,
    kTpcTof,
    nPidDetectors
  };

  std::vector<std::vector<Variable>> varPid;
  varPid.resize(nPidDetectors);
  for(auto& vP : varPid) {
    vP.resize(nProngSpecies);
  }

  std::vector<std::string> PidDetectors{"Tpc", "Tof", "TpcTof"};
  std::vector<std::string> ProngSpecies{"p", "K", "#pi"};

  for(int iDet=0; iDet<nPidDetectors; iDet++) {
    varPid.at(iDet).at(kProton) = Variable("nSig" + PidDetectors.at(iDet) + "Pr",
                                           {{"Candidates", "Lite_fNSig" + PidDetectors.at(iDet) + "Pr0"}, {"Candidates", "Lite_fNSig" + PidDetectors.at(iDet) + "Pr2"}, {"Candidates", "Lite_fCandidateSelFlag"}},
                                           []( std::vector<double>& var ) { return var.at(2)==1 ? var.at(0) : var.at(2)==2 ? var.at(1) : UndefValueFloat; });

    varPid.at(iDet).at(kPion) = Variable("nSig" + PidDetectors.at(iDet) + "Pi",
                                         {{"Candidates", "Lite_fNSig" + PidDetectors.at(iDet) + "Pi0"}, {"Candidates", "Lite_fNSig" + PidDetectors.at(iDet) + "Pi2"}, {"Candidates", "Lite_fCandidateSelFlag"}},
                                         []( std::vector<double>& var ) { return var.at(2)==1 ? var.at(1) : var.at(2)==2 ? var.at(0) : UndefValueFloat; });

    varPid.at(iDet).at(kKaon) = Variable::FromString("Candidates.Lite_fNSig" + PidDetectors.at(iDet) + "Ka1");
    varPid.at(iDet).at(kKaon).SetName("nSig" + PidDetectors.at(iDet) + "Ka");
  }

  for(auto& dt : datatypes) {
    for(auto& dcaSel : dcaFitterSelections) {
      std::vector<SimpleCut> allCuts = topoSelectionCuts;
      allCuts.emplace_back(dt);
      allCuts.emplace_back(dcaSel);

      Cuts* cutsTotal = new Cuts(dt.GetTitle() + "_" + dcaSel.GetTitle() + "_total", allCuts);

      for(int iDet=0; iDet<nPidDetectors; iDet++) {
        for(int iPs=0; iPs<nProngSpecies; iPs++) {
          task.AddH1(varPid.at(iDet).at(iPs).GetName() + "_" + dt.GetTitle() + "_" + dcaSel.GetTitle(), {"#sigma_{" + PidDetectors.at(iDet) + "} {" + ProngSpecies.at(iPs) + "}", varPid.at(iDet).at(iPs), {400, -20, 20}}, cutsTotal);
        } // nProngSpecies
      } // nPidDetectors
    } // dcaFitterSelections
  } // datatypes
}

void TopoQA(QA::Task& task) {
  struct Quantity {
    std::string name_;
    std::string name_in_tree_;
    std::string xaxis_title_;
    TAxis xaxis_;
    std::vector<SimpleCut> slice_cuts_;
  };

  auto pTCuts = HelperFunctions::CreateSliceCuts({0.f, 2.f, 4.f, 6.f, 8.f, 12.f, 16.f, 24.f}, "pT_", "Candidates.KF_fPt");

  std::vector<Quantity> vars {
    {"Mass",         "KF_fMassInv",                 "m_{pK#pi}, GeV/c^{2}", {600, 1.98, 2.58}, pTCuts},
    {"P",            "KF_fP",                       "p, GeV/c",             {600, 0,    16  }, {}    },
    {"Pt",           "KF_fPt",                      "p_{T}, GeV/c",         {600, 0,    16  }, {}    },
    {"Chi2prim_p",   "KF_fChi2PrimProton",          "#chi^{2}_{prim}{p}",   {100, 0,    60  }, {}    },
    {"Chi2prim_K",   "KF_fChi2PrimKaon",            "#chi^{2}_{prim}{K}",   {100, 0,    60  }, {}    },
    {"Chi2prim_pi",  "KF_fChi2PrimPion",            "#chi^{2}_{prim}{#pi}", {100, 0,    60  }, {}    },
    {"Chi2geo_p_pi", "KF_fChi2geoProtonPion",       "#chi^{2}_{geo}{p#pi}", {100, 0,    4   }, {}    },
    {"Chi2geo_p_K",  "KF_fChi2geoProtonKaon",       "#chi^{2}_{geo}{pK}",   {100, 0,    4   }, {}    },
    {"Chi2geo_K_pi", "KF_fChi2geoPionKaon",         "#chi^{2}_{geo}{K#pi}", {100, 0,    4   }, {}    },
    {"DCA_p_pi",     "KF_fDCAProtonPion",           "DCA{p#pi}, cm",        {100, 0,    0.1 }, {}    },
    {"DCA_p_K",      "KF_fDCAProtonKaon",           "DCA{pK}, cm",          {100, 0,    0.1 }, {}    },
    {"DCA_K_pi",     "KF_fDCAPionKaon",             "DCA{K#pi}, cm",        {100, 0,    0.1 }, {}    },
    {"Chi2geo",      "KF_fChi2geo",                 "#chi^{2}_{geo}",       {100, 0,    10  }, {}    },
    {"Chi2topo",     "KF_fChi2topo",                "#chi^{2}_{topo}",      {100, 0,    20  }, {}    },
    {"LdL",          "KF_fLdL",                     "L/#Delta L",           {100, 0,    10  }, {}    },
    {"L",            "KF_fL",                       "L, cm",                {100, 0,    0.5 }, {}    },
    {"T",            "KF_fT",                       "T, ps",                {400, 0,    5   }, {}    },
    {"isSel",        "KF_fIsSelected",              "isSel",                {10 , -5,   5   }, {}    },
    {"nPCPV",        "Lite_fNProngsContributorsPV", "nPCPV",                {6,   -1  , 5   }, {}    },
  };

  for(auto& dt : datatypes) {
    for(auto& dcaSel : dcaFitterSelections) {
      std::vector<SimpleCut> allCuts = topoSelectionCuts;
      allCuts.emplace_back(dt);
      allCuts.emplace_back(dcaSel);

      Cuts* cutsTotal = new Cuts(dt.GetTitle() + "_" + dcaSel.GetTitle() + "_total", allCuts);

      for(auto& var : vars) {
        task.AddH1(var.name_ + "_" + dt.GetTitle() + "_" + dcaSel.GetTitle(), {var.xaxis_title_, Variable::FromString("Candidates." + var.name_in_tree_), var.xaxis_}, cutsTotal);
        for(auto& slc : var.slice_cuts_) {
          std::vector<SimpleCut> allCutsSlice = allCuts;
          allCutsSlice.emplace_back(slc);
          Cuts* cutSlice = new Cuts(dt.GetTitle() + "_" + dcaSel.GetTitle() + "_" + slc.GetTitle(), allCutsSlice);

          task.AddH1(var.name_ + "_" + dt.GetTitle() + "_" + dcaSel.GetTitle() + "_" + slc.GetTitle(), {var.xaxis_title_, Variable::FromString("Candidates." + var.name_in_tree_), var.xaxis_}, cutSlice);
        } // var.slice_cuts_
      } // vars
    } // dcaFitterSelections
  } // datatypes
}

int main(int argc, char* argv[]){
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./treeKF_qa filelistname" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string filelistname = argv[1];
  treeKF_qa(filelistname);

  return 0;
}
