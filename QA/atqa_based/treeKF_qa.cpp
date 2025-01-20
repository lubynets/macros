#include <string>

#include "AnalysisTree/TaskManager.hpp"
#include "AnalysisTree/Variable.hpp"
#include "Task.hpp"

using namespace AnalysisTree;

template<typename T>
std::string to_string_with_precision(const T a_value, const int n);

std::vector<SimpleCut> CreateSliceCuts(const std::vector<float>& ranges, const std::string& cutNamePrefix, const std::string& branchFieldName);

void treeKF_qa(const std::string& filelist){
  const float HugeValue = 1e9;

  auto* man = TaskManager::GetInstance();

  auto* task = new QA::Task;
  task->SetOutputFileName("treeKF_qa.root");

  std::vector<SimpleCut> datatypes {
//     EqualsCut("Candidates.KF_fSigBgStatus", 0,    "background"),
//     EqualsCut("Candidates.KF_fSigBgStatus", 1,    "prompt"),
//     EqualsCut("Candidates.KF_fSigBgStatus", 2,    "nonprompt"),
//     EqualsCut("Candidates.KF_fSigBgStatus", 3,    "wrongswap"),
    EqualsCut("Candidates.KF_fSigBgStatus", -999, "data"),
//     SimpleCut({"Candidates.KF_fSigBgStatus"}, [](std::vector<double> par){ return par[0] != 0 && par[0] != 1 && par[0] != 2 && par[0] != 3 && par[0] != -999; }, "impossible"),
  };

  SimpleCut simpleCutSelected = RangeCut("Candidates.KF_fIsSelected", 0.9, 1.1);
//   SimpleCut simpleCutSelected = RangeCut("Candidates.KF_fIsSelected", -0.1, 1.1);

  struct Quantity {
    std::string name_;
    std::string name_in_tree_;
    std::string xaxis_title_;
    TAxis xaxis_;
    std::vector<SimpleCut> slice_cuts_;
  };

  auto pTCuts = CreateSliceCuts({0.f, 2.f, 4.f, 6.f, 8.f, 12.f, 16.f, 24.f}, "pT_", "Candidates.KF_fPt");

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

  std::vector<SimpleCut> topoSelectionCuts{
    RangeCut("Candidates.KF_fChi2PrimProton",    5,          HugeValue),
    RangeCut("Candidates.KF_fChi2PrimKaon",      5,          HugeValue),
    RangeCut("Candidates.KF_fChi2PrimPion",      5,          HugeValue),
    RangeCut("Candidates.KF_fChi2geoPionKaon",   -HugeValue, 3        ),
    RangeCut("Candidates.KF_fChi2geoProtonKaon", -HugeValue, 3        ),
    RangeCut("Candidates.KF_fChi2geoProtonPion", -HugeValue, 3        ),
    RangeCut("Candidates.KF_fDCAPionKaon",       -HugeValue, 0.01     ),
    RangeCut("Candidates.KF_fDCAProtonKaon",     -HugeValue, 0.01     ),
    RangeCut("Candidates.KF_fDCAProtonPion",     -HugeValue, 0.01     ),
    RangeCut("Candidates.KF_fChi2geo",           -HugeValue, 3        ),
    RangeCut("Candidates.KF_fChi2topo",          -HugeValue, 5        ),
    RangeCut("Candidates.KF_fLdL",               3,          HugeValue),
  };

  for(auto& dt : datatypes) {
    std::vector<SimpleCut> allCuts = topoSelectionCuts;
    allCuts.emplace_back(simpleCutSelected);
    allCuts.emplace_back(dt);

    Cuts* cutsTotal = new Cuts(dt.GetTitle() + "_total", allCuts);

    for(auto& var : vars) {
      task->AddH1(var.name_ + "_" + dt.GetTitle(), {var.xaxis_title_, Variable::FromString("Candidates." + var.name_in_tree_), var.xaxis_}, cutsTotal);
      for(auto& slc : var.slice_cuts_) {
        std::vector<SimpleCut> allCutsSlice = allCuts;
        allCutsSlice.emplace_back(slc);
        Cuts* cutSlice = new Cuts(dt.GetTitle() + "_" + slc.GetTitle(), allCutsSlice);

        task->AddH1(var.name_ + "_" + dt.GetTitle() + "_" + slc.GetTitle(), {var.xaxis_title_, Variable::FromString("Candidates." + var.name_in_tree_), var.xaxis_}, cutSlice);
      }
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
    std::cout << " ./treeKF_qa filelistname" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string filelistname = argv[1];
  treeKF_qa(filelistname);

  return 0;
}

template<typename T>
std::string to_string_with_precision(const T a_value, const int n) {
  std::ostringstream out;
  out.precision(n);
  out << std::fixed << a_value;
  return out.str();
}

std::vector<SimpleCut> CreateSliceCuts(const std::vector<float>& ranges, const std::string& cutNamePrefix, const std::string& branchFieldName) {
  std::vector<SimpleCut> sliceCuts;
  for(int iRange=0; iRange<ranges.size()-1; iRange++) {
    const std::string cutName = cutNamePrefix + to_string_with_precision(ranges.at(iRange), 2) + "_" + to_string_with_precision(ranges.at(iRange+1), 2);
    sliceCuts.emplace_back(RangeCut(branchFieldName, ranges.at(iRange), ranges.at(iRange+1), cutName));
  }

  return sliceCuts;
}
