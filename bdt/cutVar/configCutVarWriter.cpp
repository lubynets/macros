#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

template<typename T>
std::string to_string_with_precision(T a_value, int n=2);

int main(int argc, char* argv[]) {
  if(argc<2) {
    throw std::runtime_error("Not enough arguments!");
  }

  const std::string fileInName = argv[1];

  std::fstream fileIn(fileInName);

  if(!fileIn.is_open()) {
    throw std::runtime_error("Input file name is wrong!");
  }

  const std::string rawYieldName = "RawYields_Lc";
  const std::string effName = "Eff_times_Acc_Lc";

  const std::string rawYieldTemplate = "RAW_YIELD_FILES";
  const std::string effTemplate = "EFFICIENCY_FILES";
  const std::string ptBinToProcessTemplate = "PT_BIN_TO_PROCESS";

  const std::vector<std::pair<int, std::vector<double>>> fileScores {
    // data
    // gaus_poly2_fixParams
    {1, {0.03, 0.11, 0.13, 0.15, 0.17, 0.20, 0.23, 0.26, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60}}, // ct_1
    {2, {0.06, 0.18, 0.21, 0.26, 0.31, 0.34, 0.38, 0.44, 0.50, 0.57, 0.63, 0.75}}, // ct_2
    {3, {0.05, 0.23, 0.29, 0.35, 0.39, 0.45, 0.49, 0.55, 0.59, 0.65, 0.69, 0.74, 0.79}}, // ct_3
    {4, {0.10, 0.18, 0.24, 0.29, 0.33, 0.36, 0.40, 0.42, 0.47, 0.68, 0.75, 0.85}}, // ct_4
    {5, {0.05, 0.30, 0.35, 0.40, 0.45, 0.50, 0.53, 0.56, 0.60, 0.63, 0.66, 0.70, 0.80, 0.85, 0.90}}, // ct_5
  };

  std::vector<std::ofstream> fileOuts;
  for(const auto& fS : fileScores) {
    const std::string fileOutName = fileInName.substr(0, fileInName.size()-5) + "_ct" + std::to_string(fS.first) + ".json";
    fileOuts.emplace_back(fileOutName);
  } // fileScores

  std::string line;
  while (std::getline(fileIn, line)) {
    auto CopyScoreWiseLines = [&](const std::string& replaceTemplate, const std::string& fileName) {
      if(line.find(replaceTemplate) == std::string::npos) return;
      for(size_t iFO=0, nFOs=fileOuts.size(); iFO<nFOs; ++iFO) {
        for(size_t iScore=0, nScores=fileScores.at(iFO).second.size(); iScore<nScores; ++iScore) {
          const double score = fileScores.at(iFO).second.at(iScore);
          fileOuts.at(iFO) << "            \"" + fileName + ".NPgt" + to_string_with_precision(score, 2) + ".root\"";
          if(iScore != nScores-1) fileOuts.at(iFO) << ",";
          fileOuts.at(iFO) << "\n";
        } // fileScores
      } // fileOuts
    };

    CopyScoreWiseLines(rawYieldTemplate, rawYieldName);
    CopyScoreWiseLines(effTemplate, effName);

    if(line.find(ptBinToProcessTemplate) != std::string::npos) {
      for(size_t iFO=0, nFOs=fileOuts.size(); iFO<nFOs; ++iFO) {
        fileOuts.at(iFO) << "    \"pt_bin_to_process\": " << std::to_string(fileScores.at(iFO).first) << ",\n";
      } // fileOuts
    }

    if(line.find(rawYieldTemplate) == std::string::npos && line.find(effTemplate) == std::string::npos && line.find(ptBinToProcessTemplate) == std::string::npos) {
      for(auto& fO : fileOuts) {
        fO << line << "\n";
      }
    }
  }
  return 0;
}

template<typename T>
std::string to_string_with_precision(T a_value, int n) {
  std::ostringstream out;
  out.precision(n);
  out << std::fixed << a_value;
  return out.str();
}
