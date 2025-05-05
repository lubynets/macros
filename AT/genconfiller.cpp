//
// Created by oleksii on 12.04.25.
//
#include "GenericContainerFiller.hpp"

#include <iostream>

using namespace AnalysisTree;

void RunGenericContainerFiller(const std::string& fileName, int nEntries) {
  GenericContainerFiller gcf(fileName, "plainTree");
  gcf.SetNChannelsPerEntry(1);
  gcf.Run(nEntries);
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./genconfiller fileName (nEntries=ALL)" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileName = argv[1];
  const int nEntries = argc > 2 ? atoi(argv[2]) : -1;
  RunGenericContainerFiller(fileName, nEntries);

  return 0;
}