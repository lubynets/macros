#include "HelperGeneral.hpp"

#include <TH1.h>
#include <TH2.h>

#include <exception>
#include <iostream>
#include <string>
#include <vector>

using namespace HelperGeneral;

TH1* MultResponseByGenYield(const TH2* hResp, const TH1* hGen);

void relative_efficiency(const std::string& fileName) {
  TFile* fileIn = OpenFileWithNullptrCheck(fileName);

  const std::vector<std::string> promptnesses{"prompt", "nonprompt"};
  const std::vector<std::string> weightPresences{"", "_W"};

  // ========================= Configuration =================================
  const std::vector<double> pTRanges = {1, 2, 3, 4, 5, 8, 12, 20};

  std::vector<float> bdtScores{0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90};
//   for(int i=0; i<=99; i++) {
//     bdtScores.emplace_back(0.01 * i);
//   }
// //   bdtScores.emplace_back(0.01);

  const float refBdtScore{bdtScores.front()};

  if(std::find(bdtScores.begin(), bdtScores.end(), refBdtScore) == bdtScores.end()) throw std::runtime_error("relative_efficiency(): refBdtScore is not present among bdtScores");

  std::vector<std::pair<double, double>> pTIntervals{};
  for(size_t iPt = 0, nPts = pTRanges.size() - 1; iPt < nPts; ++iPt) {
//     pTIntervals.emplace_back(std::make_pair(pTRanges.at(iPt), pTRanges.at(iPt+1)));
  }

  pTIntervals.emplace_back(std::make_pair(3., 20.));

  // ==========================================================================
  auto PtRangeString = [] (const std::pair<double, double>& pTInterval) {
    return "pT_" + HelperGeneral::to_string_with_precision(pTInterval.first, 0) + "_" + HelperGeneral::to_string_with_precision(pTInterval.second, 0);
  };

  bool isFirstLoopIteration{true};
  for(const auto& promptness : promptnesses) {
    for(const auto& weightPresence : weightPresences) {
      for (size_t iPt = 0, nPts = pTIntervals.size(); iPt < nPts; ++iPt) {
        TH1* hGen = GetObjectWithNullptrCheck<TH1>(fileIn, "yields/" + promptness + "/" + PtRangeString(pTIntervals.at(iPt)) + "/gen" + weightPresence);
        TH2* hRespRef = GetObjectWithNullptrCheck<TH2>(fileIn, "effs/" + promptness + "/" + PtRangeString(pTIntervals.at(iPt)) + "/r_NPgt" + to_string_with_precision(refBdtScore, 2) + weightPresence);
        CheckHistogramsForAxisIdentity<TH2, TH2>(hRespRef, nullptr, "XY");
        CheckHistogramsForAxisIdentity(hRespRef, hGen, "X");
        TH1* hRecRef = MultResponseByGenYield(hRespRef, hGen);
        const std::string openOption = isFirstLoopIteration ? "recreate" : "update";
        for (const auto& score : bdtScores) {
          TH2* hRespScore = GetObjectWithNullptrCheck<TH2>(fileIn, "effs/" + promptness + "/" + PtRangeString(pTIntervals.at(iPt)) + "/r_NPgt" + to_string_with_precision(score, 2) + weightPresence);
          CheckHistogramsForAxisIdentity<TH2, TH2>(hRespScore, nullptr, "XY");
          CheckHistogramsForAxisIdentity(hRespScore, hGen, "X");
          TFile* fileOutScore = TFile::Open(("RelEff_Lc.NPgt" + to_string_with_precision(score, 2) + ".root").c_str(), openOption.c_str());
          TH1* hRecScore = MultResponseByGenYield(hRespScore, hGen);
          hRecScore->Divide(hRecScore, hRecRef, 1., 1., "B");
          CD(fileOutScore, PtRangeString(pTIntervals.at(iPt)));
          hRecScore->Write((promptness + weightPresence).c_str());
          fileOutScore->Close();
        } // bdtScores
        isFirstLoopIteration = false;
      } // pTRanges
    } // weightPresences
  } // promptnesses

  fileIn->Close();
}

TH1* MultResponseByGenYield(const TH2* hResp, const TH1* hGen) {
  TH1* hRec = dynamic_cast<TH1*>(hGen->Clone());
  hRec->Reset();
  hRec->GetXaxis()->SetTitle(hResp->GetXaxis()->GetTitle());
  const int nBins = hGen->GetNbinsX();

  for(int recBin=1; recBin<=nBins; ++recBin) {
    double value{};
    double err2{};
    for(int genBin=1; genBin<=nBins; ++genBin) {
      const double R = hResp->GetBinContent(recBin, genBin);
      const double T = hGen->GetBinContent(genBin);
      const double eR = hResp->GetBinError(recBin, genBin);
      const double eT = hGen->GetBinError(genBin);
      value += R * T;
      err2 += T*T*eR*eR + R*R*eT*eT;
    } // genBin
    hRec->SetBinContent(recBin, value);
    hRec->SetBinError(recBin, std::sqrt(err2));
  } // recBin

  return hRec;
}

int main(int argc, char* argv[]) {
  const int nArgs{1};
  if (argc < nArgs+1 || argc > nArgs+1) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./relative_efficiency fileName" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileName = argv[1];

  relative_efficiency(fileName);

  return 0;
}
