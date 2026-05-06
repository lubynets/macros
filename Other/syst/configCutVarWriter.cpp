#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

template<typename T>
std::string to_string_with_precision(T a_value, int n=2);

int main(int argc, char* argv[]) {
  if(argc<3) {
    throw std::runtime_error("Not enough arguments!");
  }

  const std::string fileInName = argv[1];
  const std::string rawYieldPathSuffix = argv[2];

  std::fstream fileIn(fileInName);

  if(!fileIn.is_open()) {
    throw std::runtime_error("Input file name is wrong!");
  }

  const std::string rawYieldName = "RawYields_Lc";
  const std::string effName = "Eff_times_Acc_Lc";

  const std::string rawYieldTemplate = "RAW_YIELD_FILES";
  const std::string effTemplate = "EFFICIENCY_FILES";
  const std::string ptBinToProcessTemplate = "PT_BIN_TO_PROCESS";
  const std::string rawYieldPathSuffixTemplate = "INPUT_DIR";


  const std::string rawYieldPathPrefix = "/lustre/alice/users/lubynets/syst/multiFit/rawy";
  const std::vector<std::pair<int, std::vector<double>>> fileScores {
    {1, {0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80}}, // ct_1
    {2, {0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80}}, // ct_2
    {3, {0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90}}, // ct_3
    {4, {0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90}}, // ct_4
    {5, {0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90}}, // ct_5
    {6, {0.30, 0.40, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90}}, // ct_6
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

    for(size_t iFO=0, nFOs=fileOuts.size(); iFO<nFOs; ++iFO) {
      if(line.find(ptBinToProcessTemplate) != std::string::npos) {
        fileOuts.at(iFO) << "    \"pt_bin_to_process\": " << std::to_string(fileScores.at(iFO).first) << ",\n";
      }
      if(line.find(rawYieldPathSuffixTemplate) != std::string::npos) {
        fileOuts.at(iFO) << "    \"inputdir\": \"" << rawYieldPathPrefix << "/" + rawYieldPathSuffix << "\",\n";
      }
    } // fileOuts

    if(line.find(rawYieldTemplate) == std::string::npos && line.find(effTemplate) == std::string::npos && line.find(ptBinToProcessTemplate) == std::string::npos && line.find(rawYieldPathSuffixTemplate) == std::string::npos) {
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
