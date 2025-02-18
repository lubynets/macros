#include "MultiPicture.h"

#include <iostream>
#include <array>

int main(int argc, char* argv[]) {

  MultiPicture mpic(3, 3);
  mpic.SetVerbose();

  std::string vertex = argv[1];

  std::string respull = argv[2];

  std::array<std::string, 3> component{"X", "Y", "Z"};
  std::array<std::string, 3> plot{"width", "meanband", "meanstat"};

  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
      mpic.SetPictureName(i, j, ("mc_qa2diff_" + component.at(i) + vertex + "_" + respull + "." + plot.at(j) + ".png").c_str());
    }
  }

  mpic.SetLeftMargins(0.);
  mpic.SetRightMargins(0.);
  mpic.SetTopMargins(0.);
  mpic.SetBottomMargins(0.);

  mpic.Run();

  return 0;
}
