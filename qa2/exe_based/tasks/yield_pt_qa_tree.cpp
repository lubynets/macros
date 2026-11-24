//
// Created by oleksii on 12.09.2025.
//

#include "HelperGeneral.hpp"

#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TKey.h>
#include <TTree.h>

#include <boost/program_options.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace HelperGeneral;

bool gIsDoWeight{false};
std::vector<double> gBdtSignalLowerValues{};

std::vector<double> pTRanges = {3, 4, 5, 8, 12, 20};
const std::vector<double> bdtBgUpperValuesVsPt = {0.02, 0.02, 0.02, 0.04, 0.08}; // standard
// const std::vector<double> bdtBgUpperValuesVsPt = {0.01, 0.01, 0.01, 0.02, 0.04}; // tight
// const std::vector<double> bdtBgUpperValuesVsPt = {0.03, 0.03, 0.03, 0.06, 0.12}; // loose
// const std::vector<double> bdtBgUpperValuesVsPt = {0.015, 0.015, 0.015, 0.03, 0.06}; // semitight
// const std::vector<double> bdtBgUpperValuesVsPt = {0.05, 0.05, 0.05, 0.08, 0.16}; // veryloose
// const std::vector<double> bdtBgUpperValuesVsPt = {0.025, 0.025, 0.025, 0.05, 0.10}; // semiloose
const std::vector<double> bdtNpUpperValuesVsPt = {1.00, 1.00, 1.00, 1.00, 1.00}; // set pT-dependent values

const std::string fileOutName{"yield_pt_qa_tree.root"};
const std::string promptWeightHistoName{"histoWeightcent"};
const std::string nonpromptWeightHistoName{"histoNPWeightcent"};
constexpr bool IsVerbose{true};

const std::vector<std::pair<std::string, int>> promptnesses {
  {"prompt", 1},
  {"nonprompt", 2}
};

struct Yield {
  double sum{0.};
  double sumW2{0.};
  void Fill(const double weight) {
    sum += weight;
    sumW2 += weight * weight;
  }
};

std::string GetPtCutName(size_t iPt);

std::vector<TTree*> GetTrees(TFile* fileIn, const bool isRec) {
  std::vector<TTree*> trees;
  auto AddTrees = [&](TDirectory* dir) {
    TIter next(dir->GetListOfKeys());
    while(TKey* key = dynamic_cast<TKey*>(next())) {
      TTree* tree = dir->Get<TTree>(key->GetName());
      if(!tree) continue;
      if(isRec) {
        if(!tree->GetLeaf("fPt") || !tree->GetLeaf("fOriginMcRec") || !tree->GetLeaf("ML_output_Bkg") || !tree->GetLeaf("ML_output_FD")) continue;
      } else {
        if(std::string(tree->GetName()) != "O2hfcandlcfullp") continue;
        if(!tree->GetLeaf("fPt") || !tree->GetLeaf("fOriginMcGen")) continue;
      }
      trees.emplace_back(tree);
      if(IsVerbose) std::cout << "GetTrees(): found " << dir->GetPath() << "/" << tree->GetName() << "\n";
    }
  };

  AddTrees(fileIn);
  TIter next(fileIn->GetListOfKeys());
  while(TKey* key = dynamic_cast<TKey*>(next())) {
    const std::string dirName = key->GetName();
    if(dirName.rfind("DF_", 0) != 0) continue;
    TDirectory* dir = fileIn->GetDirectory(dirName.c_str());
    if(dir) AddTrees(dir);
  }
  if(trees.empty()) throw std::runtime_error("GetTrees(): no matching trees found in " + std::string(fileIn->GetName()));
  return trees;
}

void FillYield(const std::string& fileName, const bool isRec, TH1* histoWeightPrompt, TH1* histoWeightNonPrompt) {
  if(IsVerbose) std::cout << "FillYield() started: " << (isRec ? "rec" : "gen") << "\n";
  TFile* fileIn = OpenFileWithNullptrCheck(fileName);
  const auto trees = GetTrees(fileIn, isRec);
  TFile* fileOut = TFile::Open(fileOutName.c_str(), isRec ? "recreate" : "update");
  if(!fileOut || fileOut->IsZombie()) throw std::runtime_error("FillYield(): cannot open output file");

  const size_t nPts = pTRanges.size() - 1;
  const size_t nCuts = isRec ? gBdtSignalLowerValues.size() : 1;
  const size_t nWeights = gIsDoWeight ? 2 : 1;
  std::vector<std::vector<std::vector<std::vector<Yield>>>> yields(promptnesses.size(),
    std::vector<std::vector<std::vector<Yield>>>(nPts,
      std::vector<std::vector<Yield>>(nWeights, std::vector<Yield>(nCuts))));

  for(TTree* tree : trees) {
    float fPt{};
    char fOriginMc{};
    double bdtBg{}, bdtNp{};
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("fPt", 1);
    tree->SetBranchStatus(isRec ? "fOriginMcRec" : "fOriginMcGen", 1);
    if(tree->SetBranchAddress("fPt", &fPt) < 0) throw std::runtime_error("FillYield(): cannot bind fPt");
    if(tree->SetBranchAddress(isRec ? "fOriginMcRec" : "fOriginMcGen", &fOriginMc) < 0) throw std::runtime_error("FillYield(): cannot bind origin");
    if(isRec) {
      tree->SetBranchStatus("ML_output_Bkg", 1);
      tree->SetBranchStatus("ML_output_FD", 1);
      if(tree->SetBranchAddress("ML_output_Bkg", &bdtBg) < 0) throw std::runtime_error("FillYield(): cannot bind ML_output_Bkg");
      if(tree->SetBranchAddress("ML_output_FD", &bdtNp) < 0) throw std::runtime_error("FillYield(): cannot bind ML_output_FD");
    }

    for(long long iEntry=0, nEntries=tree->GetEntries(); iEntry<nEntries; ++iEntry) {
      tree->GetEntry(iEntry);
      size_t iPt{};
      for(; iPt<nPts; ++iPt) {
        if(fPt >= pTRanges.at(iPt) && fPt < pTRanges.at(iPt+1)) break;
      }
      if(iPt == nPts) continue;
      for(size_t iPr=0; iPr<promptnesses.size(); ++iPr) {
        if(static_cast<int>(fOriginMc) != promptnesses.at(iPr).second) continue;
        if(isRec && (bdtBg >= bdtBgUpperValuesVsPt.at(iPt) || bdtNp >= bdtNpUpperValuesVsPt.at(iPt))) continue;
        const TH1* histoWeight = iPr == 0 ? histoWeightPrompt : histoWeightNonPrompt;
        // Both weights use the candidate fPt: the parent-beauty-hadron pT is absent from these trees.
        const double weight = gIsDoWeight ? histoWeight->GetBinContent(histoWeight->FindBin(fPt)) : 1.;
        if(gIsDoWeight && !std::isfinite(weight)) throw std::runtime_error("FillYield(): non-finite pT weight");
        for(size_t iB=0; iB<nCuts; ++iB) {
          if(isRec) {
            const double bsc = gBdtSignalLowerValues.at(iB);
            if(bsc >= bdtNpUpperValuesVsPt.at(iPt)) break;
            if(bdtNp <= bsc) continue;
          }
          yields.at(iPr).at(iPt).at(0).at(iB).Fill(1.);
          if(gIsDoWeight) yields.at(iPr).at(iPt).at(1).at(iB).Fill(weight);
        } // gBdtSignalLowerValues
      } // promptnesses
    } // entries
    tree->ResetBranchAddresses();
    if(IsVerbose) std::cout << "FillYield(): processed " << tree->GetName() << " (" << tree->GetEntries() << " entries)\n";
  } // trees

  for(size_t iPt=0; iPt<nPts; ++iPt) {
    for(size_t iPr=0; iPr<promptnesses.size(); ++iPr) {
      const std::string dirName = (isRec ? "rec/" : "gen/") + promptnesses.at(iPr).first + "/" + GetPtCutName(iPt);
      for(size_t iW=0; iW<nWeights; ++iW) {
        const std::string histoNameSuffix = iW == 0 ? "" : "_W";
        for(size_t iB=0; iB<nCuts; ++iB) {
          if(isRec && gBdtSignalLowerValues.at(iB) >= bdtNpUpperValuesVsPt.at(iPt)) break;
          const std::string histoName = (isRec ? "hYield_NPgt" + to_string_with_precision(gBdtSignalLowerValues.at(iB), 2) : "hYield") + histoNameSuffix;
          TH1D histoYield(histoName.c_str(), histoName.c_str(), 1, pTRanges.at(iPt), pTRanges.at(iPt+1));
          const auto& yield = yields.at(iPr).at(iPt).at(iW).at(iB);
          histoYield.SetBinContent(1, yield.sum);
          histoYield.SetBinError(1, std::sqrt(yield.sumW2));
          CD(fileOut, dirName);
          histoYield.Write();
        } // gBdtSignalLowerValues
      } // weightsPresences
    } // promptnesses
  } // pTRanges

  fileOut->Close();
  fileIn->Close();
  if(IsVerbose) std::cout << "FillYield() finished\n";
}

std::string GetPtCutName(size_t iPt) {
  return "pT_" + to_string_with_precision(pTRanges.at(iPt), 0) + "_" + to_string_with_precision(pTRanges.at(iPt+1), 0);
}

int main(int argc, char* argv[]) {
  using namespace boost::program_options;
  options_description desc ("Allowed options");
  desc.add_options()
    ("help,h", "Print usage message")
    ("i-rec", value<std::string>()->required(), "reconstructed ROOT file")
    ("i-gen", value<std::string>()->required(), "generated ROOT file")
    ("weights,w", value<std::string>(), "optional ROOT file with pT weights")
  ;

  variables_map args;
  store(parse_command_line(argc, argv, desc), args);
  if(args.count("help") || argc < 2) {
    std::cout << desc << "\n";
    return 0;
  }
  notify(args);

  if(bdtBgUpperValuesVsPt.size() != pTRanges.size()-1) throw std::runtime_error("bdtBgUpperValuesVsPt.size() != pTRanges.size()-1");
  if(bdtNpUpperValuesVsPt.size() != pTRanges.size()-1) throw std::runtime_error("bdtNpUpperValuesVsPt.size() != pTRanges.size()-1");

  for(int iB=0; iB<=99; ++iB) {
    gBdtSignalLowerValues.emplace_back(0.01 * iB);
  }

  gIsDoWeight = args.count("weights");
  TFile* fileWeight = gIsDoWeight ? OpenFileWithNullptrCheck(args["weights"].as<std::string>()) : nullptr;
  TH1* histoWeightPrompt = gIsDoWeight ? GetObjectWithNullptrCheck<TH1>(fileWeight, promptWeightHistoName) : nullptr;
  TH1* histoWeightNonPrompt = gIsDoWeight ? GetObjectWithNullptrCheck<TH1>(fileWeight, nonpromptWeightHistoName) : nullptr;

  FillYield(args["i-rec"].as<std::string>(), true, histoWeightPrompt, histoWeightNonPrompt);
  FillYield(args["i-gen"].as<std::string>(), false, histoWeightPrompt, histoWeightNonPrompt);

  if(gIsDoWeight) fileWeight->Close();
  return 0;
}
