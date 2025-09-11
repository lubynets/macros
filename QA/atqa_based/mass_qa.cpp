//
// Created by oleksii on 18.03.25.
//

#include "Task.hpp"

#include "AnalysisTree/HelperFunctions.hpp"
#include "AnalysisTree/TaskManager.hpp"

using namespace AnalysisTree;

const std::vector<float> lifetimeRanges = {0.2, 0.35, 0.5, 0.7, 0.9, 1.6};
const std::vector<float> pTRanges = {0, 2, 5, 8, 12, 20};
const TAxis massAxis = {600, 1.98, 2.58};
const std::string massAxisTitle = "m_{pK#pi} (GeV/#it{c}^{2})";
const std::pair<float, float> rapidityRanges{-0.8, 0.8};

std::string GetPtCutName(size_t iPt) {
  std::pair<size_t, size_t> iPTMinMax = (iPt == pTRanges.size()-1) ? std::pair<size_t, size_t>{0, pTRanges.size()-1} : std::pair<size_t, size_t>{iPt, iPt+1};
  return "pT_" + HelperFunctions::ToStringWithPrecision(pTRanges.at(iPTMinMax.first), 0) + "_" + HelperFunctions::ToStringWithPrecision(pTRanges.at(iPTMinMax.second), 0);
}

enum BdtClass : short {
  kBackground = 0,
  kPrompt,
  kNonPrompt,
  nBdtClasses
};

const std::vector<short> signalTargetsAtBdt {
    kPrompt,
    kNonPrompt,
};

const std::string recBranchName = "Candidates";
const std::vector<std::string> bdtClasses{"fLiteMlScoreFirstClass", "fLiteMlScoreSecondClass", "fLiteMlScoreThirdClass"};

// auto TCuts = HelperFunctions::CreateRangeCuts(lifetimeRanges, "T_", recBranchName + ".fKFT");
SimpleCut rapidityCut = RangeCut(Variable::FromString(recBranchName + ".fLiteY"), rapidityRanges.first, rapidityRanges.second);

Variable properLifetime("properLifetime", {{recBranchName, "fLiteCt"}}, [](const std::vector<double>& v) { return 100./2.997*v.at(0); });
std::vector<SimpleCut> TCuts{
  RangeCut(properLifetime, lifetimeRanges.at(0), lifetimeRanges.at(1), ("T_" + HelperFunctions::ToStringWithPrecision(lifetimeRanges.at(0), 2) + "_" + HelperFunctions::ToStringWithPrecision(lifetimeRanges.at(1), 2)).c_str()),
  RangeCut(properLifetime, lifetimeRanges.at(1), lifetimeRanges.at(2), ("T_" + HelperFunctions::ToStringWithPrecision(lifetimeRanges.at(1), 2) + "_" + HelperFunctions::ToStringWithPrecision(lifetimeRanges.at(2), 2)).c_str()),
  RangeCut(properLifetime, lifetimeRanges.at(2), lifetimeRanges.at(3), ("T_" + HelperFunctions::ToStringWithPrecision(lifetimeRanges.at(2), 2) + "_" + HelperFunctions::ToStringWithPrecision(lifetimeRanges.at(3), 2)).c_str()),
  RangeCut(properLifetime, lifetimeRanges.at(3), lifetimeRanges.at(4), ("T_" + HelperFunctions::ToStringWithPrecision(lifetimeRanges.at(3), 2) + "_" + HelperFunctions::ToStringWithPrecision(lifetimeRanges.at(4), 2)).c_str()),
  RangeCut(properLifetime, lifetimeRanges.at(4), lifetimeRanges.at(5), ("T_" + HelperFunctions::ToStringWithPrecision(lifetimeRanges.at(4), 2) + "_" + HelperFunctions::ToStringWithPrecision(lifetimeRanges.at(5), 2)).c_str()),
};

void MassQABdt(QA::Task& task);

std::vector<SimpleCut> PrepareDataTypes(const std::string& varName);
std::vector<SimpleCut> datatypes = PrepareDataTypes(recBranchName + ".fKFSigBgStatus");

const short kSignal = kNonPrompt; const std::string signalShortcut = "NP";

void mass_qa(const std::string& filelistname, bool isMc) {
  if(isMc) datatypes.pop_back();
  else     datatypes.erase(datatypes.begin(), datatypes.end()-1);

  auto* man = TaskManager::GetInstance();
  const std::string fileOutName = "mass_qa.root";

  auto* task = new QA::Task;
  task->SetOutputFileName(fileOutName);

  MassQABdt(*task);

  man->AddTask(task);
  man->Init({filelistname}, {"aTree"});
  man->SetVerbosityFrequency(100);
  man->Run();
  man->Finish();

  std::vector<std::string> pTCutNames;
  for(size_t iPt=0, nPts=pTRanges.size()-1; iPt<nPts; ++iPt) {
    pTCutNames.emplace_back(GetPtCutName(iPt));
  }
  std::vector<float> bdtSignalLowerValues;
  for (int iB = 0; iB <= 99; iB++) {
    bdtSignalLowerValues.emplace_back(0.01 * iB);
  }
  TFile* fileOut = HelperFunctions::OpenFileWithNullptrCheck(fileOutName, "update");
  for(const auto& dt : datatypes) {
    for (const auto& tCut : TCuts) {
      for (const auto& bslv : bdtSignalLowerValues) {
        std::vector<std::string> histoNames;
        for (const auto& ptcn : pTCutNames) {
          histoNames.emplace_back(dt.GetTitle() + "/" + ptcn + "/" + tCut.GetTitle() + "/hM_" + signalShortcut + "gt" + HelperFunctions::ToStringWithPrecision(bslv, 2));
        }
        TH1* histoMerged = HelperFunctions::MergeHistograms(fileOut, histoNames);
        HelperFunctions::CD(fileOut, dt.GetTitle() + "/" + GetPtCutName(pTRanges.size()-1) + "/" + tCut.GetTitle());
        histoMerged->Write(("hM_" + signalShortcut + "gt" + HelperFunctions::ToStringWithPrecision(bslv, 2)).c_str());
      } // bdtSignalLowerValues
    } // TCuts
  } // datatypes
  fileOut->Close();
}

void MassQABdt(QA::Task& task) {
  const std::vector<float> bdtBgUpperValuesVsPt = {0.02, 0.02, 0.02, 0.05, 0.08};
  if(bdtBgUpperValuesVsPt.size() != pTRanges.size() - 1) throw std::runtime_error("bdtUpperValuesVsPt.size() != pTRanges.size() - 1");

  for(size_t iPt=0, nPts=pTRanges.size()-1; iPt<nPts; ++iPt) {
    const std::string pTCutName = GetPtCutName(iPt);
    SimpleCut bdtBgScoreCut = SimpleCut({recBranchName + ".fKFPt", recBranchName + "." + bdtClasses.at(kBackground)}, [=](const std::vector<double>& par) {
      return par[0] >= pTRanges.at(iPt) && par[0] < pTRanges.at(iPt + 1) && par[1] >= 0 && par[1] < bdtBgUpperValuesVsPt.at(iPt);
    }, pTCutName);
    std::vector<float> bdtSignalLowerValues;
    for (int iB = 0; iB <= 99; iB++) {
      bdtSignalLowerValues.emplace_back(0.01 * iB);
    }
    std::vector<SimpleCut> bdtSigCuts;
    for (const auto& bslv : bdtSignalLowerValues) {
      bdtSigCuts.emplace_back(RangeCut(recBranchName + "." + bdtClasses.at(kSignal), bslv, 1, signalShortcut + "gt" + HelperFunctions::ToStringWithPrecision(bslv, 2)));
    }

    for (const auto& dt : datatypes) {
      for (const auto& tCut : TCuts) {
        task.SetTopLevelDirName(dt.GetTitle() + "/" + pTCutName + "/" + tCut.GetTitle());
        for (const auto& bsc : bdtSigCuts) {
          Cuts* cutsRec = new Cuts(dt.GetTitle() + "_" + tCut.GetTitle(), {bdtBgScoreCut, bsc, rapidityCut, tCut});
          const std::string histoName = "hM_" + bsc.GetTitle();
          task.AddH1(histoName, {massAxisTitle, Variable::FromString(recBranchName + ".fLiteM"), massAxis}, cutsRec);
        }// bdtNonPromptCuts
      } // TCuts
    } // datatypes
  } // pTRanges
} // void MassQABdt()

std::vector<SimpleCut> PrepareDataTypes(const std::string& varName) {
  std::vector<SimpleCut> result {
    RangeCut(varName, -0.1, 3.1, "all"),
//     SimpleCut({varName}, [](std::vector<double> par){ return par[0] == 0 || par[0] == 1 || par[0] == 3; }, "all_wononprompt"),
//     SimpleCut({varName}, [](std::vector<double> par){ return par[0] == 0 || par[0] == 2 || par[0] == 3; }, "all_woprompt"),
//     SimpleCut({varName}, [](std::vector<double> par){ return par[0] == 0 || par[0] == 3; }, "background"),
//     EqualsCut(varName, 1,"prompt"),
//     EqualsCut(varName, 2,"nonprompt"),
//     RangeCut(varName,  0.9, 2.1, "signal"),
//     EqualsCut(varName, 3,"wrongswap"),
    EqualsCut(varName, -999, "data"),
  };

  return result;
}

int main(int argc, char* argv[]){
  if (argc < 3) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./mass_qa filelistname mcOrData" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string filelistname = argv[1];
  const bool isMc = strcmp(argv[2], "mc") == 0 ? true : strcmp(argv[2], "data") == 0 ? false : throw std::runtime_error("mass_qa::main(): argv[2] must be either 'mc' or 'data'");

  mass_qa(filelistname, isMc);

  return 0;
}
