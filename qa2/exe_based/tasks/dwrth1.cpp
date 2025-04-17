//
// Created by oleksii on 15.04.25.
//

#include "Helper.hpp"

#include <TFile.h>

#include <iostream>
#include <string>
#include <vector>

using namespace Helper;

void DWRTH1(const std::string& fileName, const std::string& mcOrData) {
  const std::vector<std::string> vars {"nSigTpcPr", "nSigTpcPi", "nSigTpcKa", "ldl", "chi2Topo", "chi2PrimPr", "chi2PrimPi", "chi2PrimKa", "chi2Geo", "mass"};
  std::vector<std::string> dataTypes;
  if(mcOrData == "mc") {
    dataTypes.emplace_back("prompt");
    dataTypes.emplace_back("nonPrompt");
  } else if (mcOrData == "data") {
    dataTypes.emplace_back("background");
  }
  const std::vector<std::string> ptCutEdges{"0", "2", "5", "8", "12", "20"};

  TFile* fileIn = OpenFileWithNullptrCheck(fileName);

  for(const auto& var : vars) {
    std::cout << var << "\n";

    for(const auto& dt : dataTypes) {
      for(int iPt=0, nPts=ptCutEdges.size()-1; iPt<nPts; iPt++) {
        const std::string histoName = dt + "/pT_" + ptCutEdges.at(iPt) + "_" + ptCutEdges.at(iPt+1) + "/" + var;
        TH1* histo = GetObjectWithNullptrCheck<TH1>(fileIn, histoName);
        auto ranges = DetermineWorkingRangesTH1(histo);
        std::cout << ranges.first << "\t" << ranges.second << "\n";
      } // ptCutEdges
    } // dataTypes
    std::cout << "\n\n";
  } // vars


}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./dwrth1 fileName mcOrData" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileName = argv[1];
  const std::string mcOrData = argv[2];
  if(mcOrData != "mc" && mcOrData != "data") throw std::runtime_error("varCorr_qa::main(): mcOrData must be either 'mc' or 'data'");

  DWRTH1(fileName, mcOrData);

  return 0;
}