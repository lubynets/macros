#include "AnalysisTree/Constants.hpp"
#include "AnalysisTree/HelperFunctions.hpp"
#include "AnalysisTree/TaskManager.hpp"
#include "AnalysisTree/Variable.hpp"
#include "Task.hpp"

#include <string>

using namespace AnalysisTree;

const std::vector<float> pTRanges = {1, 3, 5, 8, 12, 20};
const std::vector<float> bdtBgUpperValuesVsPt = {0.02, 0.02, 0.02, 0.05, 0.08};

void PidQA(QA::Task& task, const std::string& recBranchName);

void pid_qa(const std::string& filelist, const std::string& recBranchName) {
  auto* man = TaskManager::GetInstance();

  auto* task = new QA::Task;
  task->SetOutputFileName("pid_qa.root");

  PidQA(*task, recBranchName);

  man->AddTask(task);
  man->Init({filelist}, {"aTree"});
  man->SetVerbosityFrequency(100);
  man->Run();
  man->Finish();
}

void PidQA(QA::Task& task, const std::string& recBranchName) {
  if(bdtBgUpperValuesVsPt.size() != pTRanges.size() - 1) throw std::runtime_error("bdtUpperValuesVsPt.size() != pTRanges.size() - 1");

  task.SetTopLevelDirName("PidQa");
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

  const std::vector<std::string> PidDetectors{"Tpc", "Tof", "TpcTof"};
  const std::vector<std::pair<std::string, std::string>> ProngSpecies{{"Pr", "p"}, {"Ka", "K"}, {"Pi", "#pi"}};

  std::vector<Variable> varPt;
  varPt.resize(nProngSpecies);

  varPt.at(kProton) = Variable("pT_" + ProngSpecies.at(kProton).first,
                               {{recBranchName, "fLitePtProng0"}, {recBranchName, "fLitePtProng2"}, {recBranchName, "fLiteCandidateSelFlag"}},
                               [](const std::vector<double>& var) { return var.at(2) == 1 ? var.at(0) : var.at(1); });
  varPt.at(kKaon) = Variable("pT_" + ProngSpecies.at(kKaon).first,
                             {{recBranchName, "fLitePtProng1"}},
                             [](const std::vector<double>& var) { return var.at(0); });
  varPt.at(kPion) = Variable("pT_" + ProngSpecies.at(kPion).first,
                               {{recBranchName, "fLitePtProng0"}, {recBranchName, "fLitePtProng2"}, {recBranchName, "fLiteCandidateSelFlag"}},
                               [](const std::vector<double>& var) { return var.at(2) == 1 ? var.at(1) : var.at(0); });


  SimpleCut bdtBgScoreCut = SimpleCut(std::vector<std::string>{recBranchName + ".fKFPt", recBranchName + ".fLiteMlScoreFirstClass"},
                                      [&] (const std::vector<double>& var) {
                                        bool ok{false};
                                        for(int iPt=0, nPts=pTRanges.size()-1; iPt<nPts; ++iPt) {
                                          ok |= var.at(0) >= pTRanges.at(iPt) && var.at(0) < pTRanges.at(iPt+1) && var.at(1) <= bdtBgUpperValuesVsPt.at(iPt);
                                        }
                                        return ok;
                                      });

  SimpleCut signalCut = RangeCut(recBranchName + ".fKFSigBgStatus", 0.9, 2.1);
  SimpleCut backgroundCut = SimpleCut({recBranchName + ".fKFSigBgStatus"}, [&] (const std::vector<double>& var) { return var.at(0) == 0 || var.at(0) == 3; });

  Cuts* cuts = new Cuts("cuts", {bdtBgScoreCut});
//   Cuts* cuts = new Cuts("cuts", {bdtBgScoreCut, signalCut});
//   Cuts* cuts = new Cuts("cuts", {bdtBgScoreCut, backgroundCut});

  const int nBinsPt = 100;
  const double lowPt = 0.1;
  const double hiPt = 10;
  const double hiNSigmaRange = 7.;

  std::vector<double> binEdgesPt(nBinsPt+1);
  const double logLowPt = std::log10(lowPt);
  const double logHiPt = std::log10(hiPt);
  const double logStepPt = (logHiPt - logLowPt) / nBinsPt;
  for(int iBin=0; iBin<=nBinsPt; ++iBin) {
    binEdgesPt.at(iBin) = std::pow(10., logLowPt + iBin*logStepPt);
  }

  for(int iPs=0; iPs<nProngSpecies; iPs++) {
    for(int iDet=0; iDet<nPidDetectors; iDet++) {
      const double loNSigmaRange = iDet == kTpcTof ? 0. : -hiNSigmaRange;
      task.AddH2("h2nSig" + PidDetectors.at(iDet) + ProngSpecies.at(iPs).first + "VsPtProng", {"#it{p}_{T} {" + ProngSpecies.at(iPs).second + "}", varPt.at(iPs), {static_cast<int>(binEdgesPt.size())-1, binEdgesPt.data()}},
                 {"#sigma_{" + PidDetectors.at(iDet) + "} {" + ProngSpecies.at(iPs).second + "}", Variable::FromString(recBranchName + ".fKFNSig" + PidDetectors.at(iDet) + ProngSpecies.at(iPs).first), {500, loNSigmaRange, hiNSigmaRange}}, cuts);
      task.AddH2("h2nSig" + PidDetectors.at(iDet) + ProngSpecies.at(iPs).first + "VsPtMother", {"#it{p}_{T} {pK#pi}", Variable::FromString(recBranchName + ".fKFPt"), {static_cast<int>(binEdgesPt.size())-1, binEdgesPt.data()}},
                 {"#sigma_{" + PidDetectors.at(iDet) + "} {" + ProngSpecies.at(iPs).second + "}", Variable::FromString(recBranchName + ".fKFNSig" + PidDetectors.at(iDet) + ProngSpecies.at(iPs).first), {500, loNSigmaRange, hiNSigmaRange}}, cuts);
    } //nPidDetectors
  } //nProngSpecies
}

int main(int argc, char* argv[]){
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./pid_qa filelistname" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string filelistname = argv[1];
  const std::string recBranchName = argc > 2 ? argv[2] : "Candidates";
  pid_qa(filelistname, recBranchName);

  return 0;
}
