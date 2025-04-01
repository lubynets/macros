//
// Created by oleksii on 01.04.25.
//

#include "PlainTreeFiller.hpp"
#include "TaskManager.hpp"

#include <iostream>
#include <fstream>
#include <string>

void ATPlainer(const std::string& fileName) {
  std::ofstream filelist;
  filelist.open("filelist.txt");
  filelist << fileName + "\n";
  filelist.close();

  const std::array<double, 4> sidebands{2.12, 2.20, 2.38, 2.42};
  AnalysisTree::SimpleCut invMassCut = AnalysisTree::SimpleCut({"Candidates.fKFMassInv"}, [&] (const std::vector<double>& par) { return (par[0]>sidebands.at(0) && par[0]<sidebands.at(1)) || (par[0]>sidebands.at(2) && par[0]<sidebands.at(3)); });
  AnalysisTree::Cuts* sideBandCuts = new AnalysisTree::Cuts("sideBandCuts", {invMassCut});

  AnalysisTree::SimpleCut signalCut = AnalysisTree::RangeCut("Candidates.fKFSigBgStatus", 0.9, 2.1);
  AnalysisTree::Cuts* signalCuts = new AnalysisTree::Cuts("signalCuts", {signalCut});

  auto* tree_task = new AnalysisTree::PlainTreeFiller();
  tree_task->SetOutputName("PlainTree.root", "pTree");
  std::string branchname_rec = "Candidates";
  tree_task->SetInputBranchNames({branchname_rec});
  tree_task->AddBranch(branchname_rec);
  tree_task->AddBranch(branchname_rec);
  //    tree_task->AddBranchCut(signalCuts);
  const std::vector<std::string> fields_to_preserve {
    "fKFChi2PrimProton",
    "fKFChi2PrimKaon",
    "fKFChi2PrimPion",
    "fKFChi2Geo",
    "fKFChi2Topo",
    "fKFDecayLengthNormalised",
    "fLiteNSigTpcPr",
    "fLiteNSigTpcKa",
    "fLiteNSigTpcPi",
    "fKFT",
    "fKFPt",
    "fKFMassInv",
    "fKFSigBgStatus"
  };
  tree_task->SetFieldsToPreserve(fields_to_preserve);
  tree_task->SetIsPrependLeavesWithBranchName(false);

  auto* man = AnalysisTree::TaskManager::GetInstance();
  man->AddTask(tree_task);

  man->Init({"filelist.txt"}, {"aTree"});
  man->Run(-1);// -1 = all events
  man->Finish();
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./plainer fileName" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileName = argv[1];
  ATPlainer(fileName);

  return 0;
}