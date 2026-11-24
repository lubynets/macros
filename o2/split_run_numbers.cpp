#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <vector>

int main(int argc, char* argv[]) {
  if(argc < 2 || argc > 3) {
    std::cerr << "Usage: " << argv[0] << " <input_file> [fraction]\n";
    return 1;
  }

  const std::string fileInName{argv[1]};

  double fraction{0.5};
  if(argc == 3) {
    try {
      size_t pos{};
      fraction = std::stod(argv[2], &pos);
      if(pos != std::string(argv[2]).size() || !std::isfinite(fraction) || fraction < 0. || fraction > 1.) {
        std::cerr << "Fraction must be between 0 and 1\n";
        return 1;
      }
    }
    catch(const std::exception& e) {
      std::cerr << "Invalid fraction: " << argv[2] << "\n";
      return 1;
    }
  }

  std::ifstream fileIn(fileInName);
  if(!fileIn.is_open()) {
    std::cerr << "Cannot open " << fileInName << "\n";
    return 1;
  }

  std::set<int> runNumbers{};
  std::string line{};
  while(std::getline(fileIn, line)) {
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream stream(line);
    std::string token{};
    while(stream >> token) {
      if(token.size() != 6 || !std::all_of(token.begin(), token.end(), [](unsigned char c) { return c >= '0' && c <= '9'; })) continue;
      runNumbers.insert(std::stoi(token));
    }
  }
  fileIn.close();

  std::vector<int> runs(runNumbers.begin(), runNumbers.end());
  std::random_device rd{};
  std::mt19937 generator(rd());
  std::shuffle(runs.begin(), runs.end(), generator);

  const size_t nFirst = std::round(fraction * runs.size());
  std::set<int> firstSet(runs.begin(), runs.begin() + nFirst);
  std::set<int> secondSet(runs.begin() + nFirst, runs.end());

  std::ofstream fileOut1("runnumbers_1.txt");
  std::ofstream fileOut2("runnumbers_2.txt");
  if(!fileOut1.is_open() || !fileOut2.is_open()) {
    std::cerr << "Cannot create output files\n";
    return 1;
  }

  for(const auto& run : firstSet) fileOut1 << run << "\n";
  for(const auto& run : secondSet) fileOut2 << run << "\n";

  std::cout << "Total unique runs: " << runs.size() << "\n";
  std::cout << "First file: " << firstSet.size() << " runs\n";
  std::cout << "Second file: " << secondSet.size() << " runs\n";

  return 0;
}
