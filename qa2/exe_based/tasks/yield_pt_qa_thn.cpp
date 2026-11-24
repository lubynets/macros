//
// Created by oleksii on 20.11.2026.
//

#include "HelperGeneral.hpp"
#include "HelperMath.hpp"

#include <TAxis.h>
#include <THnSparse.h>
#include <TH1.h>

#include <boost/program_options.hpp>

#include <iostream>
#include <string>
#include <string_view>
#include <vector>
#include <array>
#include <algorithm>
#include <stdexcept>

using namespace HelperGeneral;
using namespace HelperMath;
using namespace std::string_literals;

std::string gCmdArgs{};
bool gIsSaveMeta{false};
bool gIsDoWeight{false};

std::vector<double> pTRanges = {1., 2., 3., 4., 6., 8., 12., 20.};
const std::vector<double> bdtBgUpperValuesVsPt = {0.02, 0.07, 0.03, 0.03, 0.03, 0.03, 0.03}; // standard
// const std::vector<double> bdtBgUpperValuesVsPt = {0.01, 0.01, 0.01, 0.02, 0.04}; // tight
// const std::vector<double> bdtBgUpperValuesVsPt = {0.03, 0.03, 0.03, 0.06, 0.12}; // loose
// const std::vector<double> bdtBgUpperValuesVsPt = {0.015, 0.015, 0.015, 0.03, 0.06}; // semitight
// const std::vector<double> bdtBgUpperValuesVsPt = {0.05, 0.05, 0.05, 0.08, 0.16}; // veryloose
// const std::vector<double> bdtBgUpperValuesVsPt = {0.025, 0.025, 0.025, 0.05, 0.10}; // semiloose

const std::vector<double> bdtNpUpperValuesVsPt = {0.80, 0.80, 0.80, 0.80, 0.80, 0.80, 0.80}; // set pT-dependent upper cuts here

constexpr std::string_view PtAxisTitle = "#it{p}_{T}(#Lambda_{c}^{+}) (GeV/#it{c})";
constexpr std::string_view PtBAxisTitle = "#it{p}_{T}^{B} (GeV/#it{c})";
constexpr std::string_view BgAxisTitle = "BDT bkg score (Lc)";
constexpr std::string_view NpAxisTitle = "BDT non-prompt score (Lc)";
constexpr std::string_view SignalTypeAxisTitle = "candidates type";

const std::string fileOutName{"yield_pt_qa_thn.root"};
constexpr bool IsVerbose{true};

const std::vector<std::pair<std::string, double>> promptnesses {
  {"prompt", 1.},
  {"nonprompt", 2.}
};

std::string GetPtCutName(size_t iPt) {
  return "pT_" + to_string_with_precision(pTRanges.at(iPt), 0) + "_" + to_string_with_precision(pTRanges.at(iPt+1), 0);
}

void FillYield(const std::string& fileName, const std::string& filePtWeightName, const bool isRec, const std::string& pCMM, const std::string& npCMM) {
  if(IsVerbose) std::cout << "FillYield() started\n";
  TFile* fileIn = OpenFileWithNullptrCheck(fileName);
  const std::string fileOpenOption = isRec ? "recreate" : "update";
  TFile* fileOut = TFile::Open(fileOutName.c_str(), fileOpenOption.c_str());
  TFile* fileWeight = gIsDoWeight ? OpenFileWithNullptrCheck(filePtWeightName) : nullptr;
  TH1* histoWeightPrompt = gIsDoWeight ? GetObjectWithNullptrCheck<TH1>(fileWeight, "histoWeight" + pCMM) : nullptr;
  TH1* histoWeightNonPrompt = gIsDoWeight ? GetObjectWithNullptrCheck<TH1>(fileWeight, "histoNPWeight" + npCMM) : nullptr;
  THnSparse* histoRecOrGen = GetObjectWithNullptrCheck<THnSparse>(fileIn, "hf-task-lc/"s + (isRec ? "hnLcVarsWithBdt" : "hnLcVarsGen"));
  const std::map<std::string_view, int> axesIndices = MapTHnSparseAxesIndices(histoRecOrGen);
  THnSparse* histoPromptWeighted = gIsDoWeight ? dynamic_cast<THnSparse*>(histoRecOrGen->Clone()) : nullptr;
  THnSparse* histoNonPromptWeighted = gIsDoWeight ? dynamic_cast<THnSparse*>(histoRecOrGen->Clone()) : nullptr;
  if (gIsDoWeight) {
    ScaleTHnSparseWithWeight(histoPromptWeighted, axesIndices.at(PtAxisTitle), histoWeightPrompt);
    ScaleTHnSparseWithWeight(histoNonPromptWeighted, axesIndices.at(PtBAxisTitle), histoWeightNonPrompt);
  }

  std::vector<double> bdtSignalLowerValues{};
  if (isRec) {
    const double bdtNpUpperMax = *std::max_element(bdtNpUpperValuesVsPt.begin(), bdtNpUpperValuesVsPt.end());
    for (int iB = 0; iB <= 99 && 0.01 * iB < bdtNpUpperMax; iB++) {
      bdtSignalLowerValues.emplace_back(0.01 * iB);
    }
    //   bdtSignalLowerValues = {0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90};
  } else {
    bdtSignalLowerValues.emplace_back(UndefValueDouble);
  }

  auto ProcessTHnSparse = [&](THnSparse* histoIn, const std::string& histoNameSuffix="", const std::vector<std::pair<std::string, double>>& promptnessesToProcess=promptnesses) {
    if(IsVerbose) std::cout << "ProcessTHnSparse() started\n";
    CheckTAxisForRanges(*histoIn->GetAxis(axesIndices.at(PtAxisTitle)), pTRanges);
    CheckTAxisForRanges(*histoIn->GetAxis(axesIndices.at(SignalTypeAxisTitle)), {1., 2., 3.});
    if(isRec) {
      CheckTAxisForRanges(*histoIn->GetAxis(axesIndices.at(BgAxisTitle)), bdtBgUpperValuesVsPt);
      CheckTAxisForRanges(*histoIn->GetAxis(axesIndices.at(NpAxisTitle)), bdtNpUpperValuesVsPt);
      CheckTAxisForRanges(*histoIn->GetAxis(axesIndices.at(NpAxisTitle)), bdtSignalLowerValues);
    }

    for(const auto& promptness : promptnessesToProcess) {
      if(IsVerbose) std::cout << "ProcessTHnSparse(): promptness = " << promptness.first << "\n";
      const std::string dirName = (isRec ? "rec/" : "gen/") + promptness.first;
      SetTHnSparseAxisRanges(histoIn, axesIndices.at(SignalTypeAxisTitle), promptness.second, promptness.second+1.f);
      // for rec - real bdtSignalLowerValues; for gen - fake 1-element vector for universality reasons
      for(const auto& bsc : bdtSignalLowerValues) {
        if(IsVerbose) std::cout << "ProcessTHnSparse(): bsc = " << bsc << "\n";
        const std::string histoName = isRec ?
                                      "hPt_NPgt" + to_string_with_precision(bsc, 2) + histoNameSuffix :
                                      "hPt" + histoNameSuffix;
        TH1D* histoYield = new TH1D(histoName.c_str(), PtAxisTitle.data(), pTRanges.size()-1, pTRanges.data());
        histoYield->SetDirectory(nullptr);

        for(size_t iPt=0, nPts=pTRanges.size()-1; iPt<nPts; ++iPt) {
          const double bdtNpUpper = bdtNpUpperValuesVsPt.at(iPt);
          if(isRec && bdtNpUpper <= bsc) continue;
          SetTHnSparseAxisRanges(histoIn, axesIndices.at(PtAxisTitle), pTRanges.at(iPt), pTRanges.at(iPt + 1));
          if(isRec) {
            SetTHnSparseAxisRanges(histoIn, axesIndices.at(BgAxisTitle), 0., bdtBgUpperValuesVsPt.at(iPt));
            SetTHnSparseAxisRanges(histoIn, axesIndices.at(NpAxisTitle), bsc, bdtNpUpper);
          }

          TH1* histoProjection = histoIn->Projection(axesIndices.at(PtAxisTitle), histoNameSuffix.empty() ? "A" : "EA");
          histoProjection->SetDirectory(nullptr);
          TH1* histoRebinned = histoProjection->Rebin(pTRanges.size()-1, (histoName + "_tmp").c_str(), pTRanges.data());
          histoRebinned->SetDirectory(nullptr);
          if(histoNameSuffix.empty()) histoRebinned->Sumw2(false);
          histoYield->Add(histoRebinned);
          delete histoRebinned;
          delete histoProjection;

          SetTHnSparseAxisRanges(histoIn, axesIndices.at(PtAxisTitle));
          if(isRec) {
            SetTHnSparseAxisRanges(histoIn, axesIndices.at(BgAxisTitle));
            SetTHnSparseAxisRanges(histoIn, axesIndices.at(NpAxisTitle));
          }
        } // pTRanges
        CD(fileOut, dirName);
        histoYield->Write(histoName.c_str());
        delete histoYield;
      } // bdtSignalLowerValues
      SetTHnSparseAxisRanges(histoIn, axesIndices.at(SignalTypeAxisTitle));
    } // promptnessesToProcess
    if(IsVerbose) std::cout << "ProcessTHnSparse() finished\n";
  };

  ProcessTHnSparse(histoRecOrGen);
  if (gIsDoWeight) {
    ProcessTHnSparse(histoPromptWeighted, "_W", {promptnesses.at(0)});
    ProcessTHnSparse(histoNonPromptWeighted, "_W", {promptnesses.at(1)});
  }

  fileOut->Close();
  fileIn->Close();
  if (gIsDoWeight) fileWeight->Close();
  delete histoPromptWeighted;
  delete histoNonPromptWeighted;
  if(IsVerbose) std::cout << "FillYield() finished\n";
}

int main(int argc, char* argv[]) {
  using namespace boost::program_options;
  options_description desc ("Allowed options");
  desc.add_options()
    ("help,h", "Print usage message")
    ("input-file,i", value<std::string>()->required(), "file with O2Physics::taskLc output (AnalysisResults.root)")
    ("weight-file,w", value<std::string>()->default_value(""), "file with pT weights")
    ("prompt-weight-suffix", value<std::string>()->default_value("cent"), "prompt weight histogram suffix")
    ("nonprompt-weight-suffix", value<std::string>()->default_value("cent"), "nonprompt weight histogram suffix")
    ("save-meta,m", value<bool>()->implicit_value(true)->default_value(false), "save meta data")
  ;

  variables_map args;
  store(parse_command_line(argc, argv, desc), args);

  if (args.count("help") || argc < 2) {
    std::cout << desc << "\n";
    return 0;
  }
  notify (args);

  for(int iArg=1; iArg<argc; ++iArg) {
    gCmdArgs += argv[iArg];
    gCmdArgs += " ";
  }

  if(bdtBgUpperValuesVsPt.size() != pTRanges.size() - 1) throw std::runtime_error("bdtBgUpperValuesVsPt.size() != pTRanges.size() - 1");
  if(bdtNpUpperValuesVsPt.size() != pTRanges.size() - 1) throw std::runtime_error("bdtNpUpperValuesVsPt.size() != pTRanges.size() - 1");

  const std::string fileNameIn = args["input-file"].as<std::string>();
  const std::string filePtWeightName = args["weight-file"].as<std::string>();
  const std::string pCMM = args["prompt-weight-suffix"].as<std::string>();
  const std::string npCMM = args["nonprompt-weight-suffix"].as<std::string>();
  gIsSaveMeta = args["save-meta"].as<bool>();

  gIsDoWeight = !filePtWeightName.empty();
  const std::string& fileName = ReadNthLine(fileNameIn);

  FillYield(fileName, filePtWeightName, true, pCMM, npCMM);
  FillYield(fileName, filePtWeightName, false, pCMM, npCMM);

  if (gIsSaveMeta) {
    TFile* fileOut = OpenFileWithNullptrCheck(fileOutName, "update");
    HelperGeneral::EmbedMetaData(fileOut, "yield_pt_qa_thn.cpp", gCmdArgs);
    fileOut->Close();
  }

  return 0;
}
