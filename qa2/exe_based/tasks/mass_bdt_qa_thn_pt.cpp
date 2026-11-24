//
// Created by oleksii on 13.11.2026.
//

#include "HelperGeneral.hpp"

#include <TAxis.h>
#include <THnSparse.h>
#include <TH1.h>

#include <boost/program_options.hpp>

#include <array>
#include <iostream>
#include <string>
#include <string_view>
#include <utility>

using namespace HelperGeneral;

constexpr std::array pTRanges = {1., 2., 3., 4., 6., 8., 12., 20.};
const std::string pTAxisTitle = "#it{p}_{T}(#Lambda_{c}^{+}) (GeV/#it{c})";

std::string gCmdArgs{};
bool gIsSaveMeta{false};

constexpr std::array bdtBgUpperValuesVsPt = {0.02, 0.07, 0.03, 0.03, 0.03, 0.03, 0.03}; // standard
// const std::vector<double> bdtBgUpperValuesVsPt = {0.01, 0.01, 0.01, 0.02, 0.04}; // tight
// const std::vector<double> bdtBgUpperValuesVsPt = {0.03, 0.03, 0.03, 0.06, 0.12}; // loose
// const std::vector<double> bdtBgUpperValuesVsPt = {0.015, 0.015, 0.015, 0.03, 0.06}; // semitight
// const std::vector<double> bdtBgUpperValuesVsPt = {0.05, 0.05, 0.05, 0.08, 0.16}; // veryloose
// const std::vector<double> bdtBgUpperValuesVsPt = {0.025, 0.025, 0.025, 0.05, 0.10}; // semiloose

constexpr std::array bdtNpUpperValuesVsPt = {0.80, 0.80, 0.80, 0.80, 0.80, 0.80, 0.80}; // set pT-dependent upper cuts here

const std::string bgAxisTitle = "BDT bkg score (Lc)";
const std::string npAxisTitle = "BDT non-prompt score (Lc)";
const std::string massAxisTitle = "inv. mass (p K #pi) (GeV/#it{c}^{2})";

constexpr int Verobsity{3};

std::string GetPtCutName(size_t iPt);

void MassBdtQaThn(const std::string& fileNameIn) {
//  LoadMacro("styles/mc_qa2.style.cc");
  const std::string fileName = ReadNthLine(fileNameIn);

  TFile* fileIn = OpenFileWithNullptrCheck(fileName);
  TFile* fileOut = TFile::Open("mass_bdt_qa_thn_pt.root", "recreate");

  THnSparse* histoIn = GetObjectWithNullptrCheck<THnSparse>(fileIn, "hf-task-lc/hnLcVarsWithBdt");

  const std::map<std::string_view, int> axesIndices = MapTHnSparseAxesIndices(histoIn);

  static_assert(bdtBgUpperValuesVsPt.size() == pTRanges.size() - 1, "bdtBgUpperValuesVsPt.size() != pTRanges.size() - 1");
  static_assert(bdtNpUpperValuesVsPt.size() == pTRanges.size() - 1, "bdtNpUpperValuesVsPt.size() != pTRanges.size() - 1");

  std::vector<double> bdtScanValues{};
  for (int iB = 0; iB <= 99; iB++) {
    bdtScanValues.emplace_back(0.01 * iB);
  }

  CheckTAxisForRanges(*histoIn->GetAxis(axesIndices.at(pTAxisTitle)), std::vector<double>(pTRanges.begin(), pTRanges.end())); // TODO vector -> array in CheckTAxisForRanges or templatize, or std::span (need c++20)
  CheckTAxisForRanges(*histoIn->GetAxis(axesIndices.at(bgAxisTitle)), std::vector<double>(bdtBgUpperValuesVsPt.begin(), bdtBgUpperValuesVsPt.end()));
  CheckTAxisForRanges(*histoIn->GetAxis(axesIndices.at(npAxisTitle)), bdtScanValues);
  CheckTAxisForRanges(*histoIn->GetAxis(axesIndices.at(npAxisTitle)), std::vector<double>(bdtNpUpperValuesVsPt.begin(), bdtNpUpperValuesVsPt.end()));

  std::vector<std::string> pTCutNames;
  for(size_t iPt=0, nPts=pTRanges.size()-1; iPt<nPts; ++iPt) {
    pTCutNames.emplace_back(GetPtCutName(iPt));
  }

  for(size_t iPt=0, nPts=pTRanges.size()-1; iPt<nPts; ++iPt) {
    if(Verobsity >= 1) std::cout << "\nProcessing iPt = " << iPt << "\n";
    SetTHnSparseAxisRanges(histoIn, axesIndices.at(pTAxisTitle), pTRanges.at(iPt), pTRanges.at(iPt + 1));
    SetTHnSparseAxisRanges(histoIn, axesIndices.at(bgAxisTitle), 0., bdtBgUpperValuesVsPt.at(iPt));
    const double bdtNpUpper = bdtNpUpperValuesVsPt.at(iPt);
    const std::string dirName = pTCutNames.at(iPt);
    if(Verobsity >= 3) std::cout << "Processing bdtScan = ";
    for (const auto& bdtScan: bdtScanValues) {
      if(bdtNpUpper <= bdtScan) continue;
      if(Verobsity >= 3) std::cout << bdtScan << " ";
      SetTHnSparseAxisRanges(histoIn, axesIndices.at(npAxisTitle), bdtScan, bdtNpUpper);
      TH1D* histoMass = histoIn->Projection(axesIndices.at(massAxisTitle));
      histoMass->SetDirectory(nullptr);
      std::string histoName = "hM_NPgt" + to_string_with_precision(bdtScan, 2) + "_NPlt" + to_string_with_precision(bdtNpUpper, 2);
      CD(fileOut, dirName);
      histoMass->Write(histoName.c_str());
      delete histoMass;
      SetTHnSparseAxisRanges(histoIn, axesIndices.at(npAxisTitle));
    } // bdtScanValues
    if(Verobsity >= 3) std::cout << "\n";
    SetTHnSparseAxisRanges(histoIn, axesIndices.at(pTAxisTitle));
    SetTHnSparseAxisRanges(histoIn, axesIndices.at(bgAxisTitle));
  } // pTRanges

  if (gIsSaveMeta) HelperGeneral::EmbedMetaData(fileOut, "mass_bdt_qa_thn_pt.cpp", gCmdArgs);
  fileOut->Close();
  fileIn->Close();
}

std::string GetPtCutName(size_t iPt) {
  return "pT_" + to_string_with_precision(pTRanges.at(iPt), 0) + "_" + to_string_with_precision(pTRanges.at(iPt+1), 0);
}

int main(int argc, char* argv[]) {
  using namespace boost::program_options;
  options_description desc ("Allowed options");
  desc.add_options()
    ("help,h", "Print usage message")
    ("input-file,i", value<std::string>()->required(), "file with O2Physics::taskLc output (AnalysisResults.root)")
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

  const std::string fileNameIn = args["input-file"].as<std::string>();
  gIsSaveMeta = args["save-meta"].as<bool>();

  MassBdtQaThn(fileNameIn);

  return 0;
}
