#include "MultiPicture.h"

int main(int argc, char* argv[]) {
  constexpr int NPts = 8;
  constexpr int NCts = 1;

  constexpr bool IsCtHorizontal{false};

  constexpr int RebinFactor = 2;

  //===========================================

  constexpr int Nx = IsCtHorizontal ? NCts : NPts;
  constexpr int Ny = IsCtHorizontal ? NPts : NCts;

  static_assert(Nx % RebinFactor == 0, "Nx % RebinFactor != 0");

  const std::string ccPrefix = argv[1];

  MultiPicture mpic(Nx/RebinFactor, Ny*RebinFactor);
  for(int iPt=0; iPt<NPts; ++iPt) {
    for(int iCt=0; iCt<NCts; ++iCt) {
      const int iX = (IsCtHorizontal ? iCt : iPt) % (Nx/RebinFactor);
      const int iY = (IsCtHorizontal ? iPt : iCt) + (IsCtHorizontal ? iCt : iPt) / (Nx/RebinFactor);
      mpic.SetPictureName(iX, iY, "process/" + ccPrefix + "_pt" + std::to_string(iPt) + "_ct" + std::to_string(iCt) + ".png");
    }
  }

  mpic.SetVerbose();
//   mpic.SetSaveIntermediatePictures();
  mpic.Run();

  return 0;
}
