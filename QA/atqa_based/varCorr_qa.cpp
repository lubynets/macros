//
// Created by oleksii on 04.04.25.
//
#include "Task.hpp"

#include "AnalysisTree/HelperFunctions.hpp"
#include "AnalysisTree/TaskManager.hpp"

#include <string>

using namespace AnalysisTree;

void VarCorrQA(QA::Task& task, const std::string& mcOrData);

void varCorr_qa(const std::string& filelist, const std::string& mcOrData, int nEntries) {
  auto* man = TaskManager::GetInstance();

  auto* task = new QA::Task;
  task->SetOutputFileName("varCorr_qa.root");

  VarCorrQA(*task, mcOrData);

  man->AddTask(task);
  man->Init({filelist}, {"aTree"});
  man->SetVerbosityPeriod(10000);
  man->Run(nEntries);
  man->Finish();
}

void VarCorrQA(QA::Task& task, const std::string& mcOrData) {
  const std::array<double, 4> sidebands{2.12, 2.20, 2.38, 2.42};
  std::vector<SimpleCut> dataTypes;
  if(mcOrData == "mc") {
    dataTypes.emplace_back(EqualsCut("Candidates.fKFSigBgStatus", 1, "prompt"));
    dataTypes.emplace_back(EqualsCut("Candidates.fKFSigBgStatus", 2, "nonPrompt"));
  } else if (mcOrData == "data") {
    dataTypes.emplace_back(SimpleCut({"Candidates.fKFMassInv"}, [=] (const std::vector<double>& par) { return (par[0]>sidebands.at(0) && par[0]<sidebands.at(1)) || (par[0]>sidebands.at(2) && par[0]<sidebands.at(3)); }, "background"));
  }

  auto pTCuts = HelperFunctions::CreateRangeCuts({0.f, 2.f, 5.f, 8.f, 12.f, 20.f}, "pT_", "Candidates.fKFPt");

  struct Quantity {
    std::string name_;
    std::string name_in_tree_;
    std::string title_;
    std::string unit_;
    TAxis axis_;
  };

  const std::vector<Quantity> vars {
    {"nSigTpcPr", "fLiteNSigTpcPr", "N#sigma_{TPC} [p]", "", {100, -5, 5}},
    {"nSigTpcKa", "fLiteNSigTpcKa", "N#sigma_{TPC} [K]", "", {100, -5, 5}},
    {"nSigTpcPi", "fLiteNSigTpcPi", "N#sigma_{TPC} [#pi]", "", {100, -5, 5}},
    {"chi2PrimPr", "fKFChi2PrimProton", "#chi^{2}_{prim} [p]", "", {225, -10, 100}},
    {"chi2PrimKa", "fKFChi2PrimKaon", "#chi^{2}_{prim} [K]", "", {225, -10, 100}},
    {"chi2PrimPi", "fKFChi2PrimPion", "#chi^{2}_{prim} [#pi]", "", {225, -10, 100}},
    {"chi2Geo", "fKFChi2Geo", "#chi^{2}_{geo} [pK#pi]", "", {225, -10, 100}},
    {"chi2Topo", "fKFChi2Topo", "#chi^{2}_{topo} [#Lambda_{c}]", "", {225, -10, 100}},
    {"ldl", "fKFDecayLengthNormalised", "L/#Delta L", "", {225, -10, 100}},
    {"mass", "fKFMassInv", "m_{pK#pi}", "GeV/#it{c}^{2}", {300, 2.12, 2.42}}
  };

  for(const auto& dt : dataTypes) {
    for(const auto& ptCut : pTCuts) {
      const std::string cutName = dt.GetTitle() + "/" + ptCut.GetTitle();
      task.SetTopLevelDirName(cutName);
      std::cout << "cutName = " << cutName << "\n";
      Cuts* cut2D = new Cuts(cutName, {dt, ptCut});
        for (int iVar = 0, nVars = vars.size(); iVar < nVars; iVar++) {
          const Quantity& xVar = vars.at(iVar);
          const std::string xVarAxisTitle = xVar.unit_.empty() ? xVar.title_ : xVar.title_ + " (" + xVar.unit_ + ")";
          task.AddH1(xVar.name_, {xVarAxisTitle, Variable::FromString("Candidates." + xVar.name_in_tree_), xVar.axis_}, cut2D);
          for (int jVar = iVar + 1; jVar < nVars; jVar++) {
            const Quantity& yVar = vars.at(jVar);
            const std::string yVarAxisTitle = yVar.unit_.empty() ? yVar.title_ : yVar.title_ + " (" + yVar.unit_ + ")";
            const std::string histoName = xVar.name_ + "_vs_" + yVar.name_;
            task.AddH2(histoName, {xVarAxisTitle, Variable::FromString("Candidates." + xVar.name_in_tree_), xVar.axis_},
                       {yVarAxisTitle, Variable::FromString("Candidates." + yVar.name_in_tree_), yVar.axis_}, cut2D);
          } // jVar : nVars
      } // iVar : nVars
    } // pTCuts
  } // dataTypes
}

int main(int argc, char* argv[]){
  if (argc < 3) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./varCorr_qa filelistname mcOrData (nEntries=ALL)" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string filelistname = argv[1];
  const std::string mcOrData = argv[2];
  if(mcOrData != "mc" && mcOrData != "data") throw std::runtime_error("varCorr_qa::main(): mcOrData must be either 'mc' or 'data'");
  const int nEntries = argc>3 ? atoi(argv[3]) : -1;
  varCorr_qa(filelistname, mcOrData, nEntries);

  return 0;
}