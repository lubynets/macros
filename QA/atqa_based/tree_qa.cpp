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

//   SimpleCut simpleCutSelected = RangeCut("Candidates.KF_fIsSelected", 0.9, 1.1);
  SimpleCut simpleCutSelected = RangeCut("Candidates.KF_fIsSelected", -0.1, 1.1);

  struct Quantity {
    std::string name_;
    std::string name_in_tree_;
    std::string xaxis_title_;
    TAxis xaxis_;
  };

  std::vector<Quantity> vars {
    {"Mass",         "KF_fMassInv",                 "m_{pK#pi}, GeV/c^{2}", {600, 1.98, 2.58}},
    {"Pt",           "KF_fPt",                      "p_{T}, GeV/c",         {600, 0,    12  }},
    {"Chi2prim_p",   "KF_fChi2PrimProton",          "#chi^{2}_{prim}{p}",   {100, 0,    60  }},
    {"Chi2prim_K",   "KF_fChi2PrimKaon",            "#chi^{2}_{prim}{K}",   {100, 0,    60  }},
    {"Chi2prim_pi",  "KF_fChi2PrimPion",            "#chi^{2}_{prim}{#pi}", {100, 0,    60  }},
    {"Chi2geo_p_pi", "KF_fChi2geoProtonPion",       "#chi^{2}_{geo}{p#pi}", {100, 0,    4   }},
    {"Chi2geo_p_K",  "KF_fChi2geoProtonKaon",       "#chi^{2}_{geo}{pK}",   {100, 0,    4   }},
    {"Chi2geo_K_pi", "KF_fChi2geoPionKaon",         "#chi^{2}_{geo}{K#pi}", {100, 0,    4   }},
    {"DCA_p_pi",     "KF_fDCAProtonPion",           "DCA{p#pi}, cm",        {100, 0,    0.1 }},
    {"DCA_p_K",      "KF_fDCAProtonKaon",           "DCA{pK}, cm",          {100, 0,    0.1 }},
    {"DCA_K_pi",     "KF_fDCAPionKaon",             "DCA{K#pi}, cm",        {100, 0,    0.1 }},
    {"Chi2geo",      "KF_fChi2geo",                 "#chi^{2}_{geo}",       {100, 0,    10  }},
    {"Chi2topo",     "KF_fChi2topo",                "#chi^{2}_{topo}",      {100, 0,    20  }},
    {"LdL",          "KF_fLdL",                     "L/#Delta L",           {100, 0,    10  }},
    {"L",            "KF_fL",                       "L, cm",                {100, 0,    0.5 }},
    {"T",            "KF_fT",                       "T, ps",                {400, 0,    5   }},
    {"nPCPV",        "Lite_fNProngsContributorsPV", "nPCPV",                {6,   -1  , 5   }},
  };

  for(auto& dt : datatypes) {
    Cuts* cut = new Cuts(dt.name_, {simpleCutSelected, dt.cut_});
    for(auto& var : vars) {
      task->AddH1((var.name_ + "_" + dt.name_).c_str(), {var.xaxis_title_, Variable::FromString(("Candidates." + var.name_in_tree_).c_str()), var.xaxis_}, cut);
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
    std::cout << " ./tree_qa filelistname" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string filelistname = argv[1];
  tree_qa(filelistname);

  return 0;
}
