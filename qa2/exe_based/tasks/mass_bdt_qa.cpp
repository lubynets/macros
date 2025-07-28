//
// Created by oleksii on 12.05.25.
//
#include "HelperGeneral.hpp"
#include "HelperMath.hpp"
#include "HelperPlot.hpp"

#include <TLegend.h>

#include <iostream>

using namespace HelperGeneral;
using namespace HelperMath;
using namespace HelperPlot;

void MassBdtQa(const std::string& fileName, const std::string& dataType) {
  LoadMacro("styles/mc_qa2.style.cc");
  TFile* fileIn = OpenFileWithNullptrCheck(fileName);

  const std::vector<Color_t> palette{kBlue, kRed, kGreen+2, kMagenta, kBlack, kOrange+2};

  const int rebinFactor = 4;

  std::vector<float> pTRanges{0.f, 2.f, 5.f, 8.f, 12.f, 20.f};
  std::vector<float> lifetimeRanges{0.2, 0.35, 0.5, 0.7, 0.9, 1.6};

  const std::string promptBdtCut = HelperGeneral::to_string_with_precision(0., 1);

  std::string priBra;
  auto PrintCanvases = [&] (const std::vector<float>& sliceRanges, const std::string& sliceVarName, const int precision) {
    priBra = "(";
    const std::string fileOutName = sliceVarName + "_Pgt" + promptBdtCut;
    for (int iSlice = 0; iSlice < sliceRanges.size() - 1; iSlice++) {
      const std::string sliceVarLow = HelperGeneral::to_string_with_precision(sliceRanges.at(iSlice), precision);
      const std::string sliceVarUp = HelperGeneral::to_string_with_precision(sliceRanges.at(iSlice + 1), precision);
      const std::string pTName = sliceVarName + "_" + sliceVarLow + "_" + sliceVarUp;
      TCanvas cc("cc", "");
      cc.SetCanvasSize(1200, 800);
      TLegend* leg = new TLegend(0.8, 0.4, 0.95, 0.8);
      leg->SetHeader(("bdt(prompt) > " + promptBdtCut).c_str());
      for (int iBg = 1; iBg <= 10; iBg++) {
        const std::string sBgCut = HelperGeneral::to_string_with_precision((11 - iBg) * 0.01, 2);
        const std::string histoName =
                dataType + "/BGlt" + sBgCut + "/" + pTName + "/hMass_BGlt" + sBgCut + "_Pgt" + promptBdtCut;
        TH1* h = GetObjectWithNullptrCheck<TH1>(fileIn, histoName);
        h->UseCurrentStyle();
        if (rebinFactor != 1) h->Rebin(rebinFactor);
        h->SetLineColor(palette.at((iBg - 1) % palette.size()));
        leg->AddEntry(h, ("bdt(BG) < " + sBgCut).c_str(), "L");
        if (iBg == 1) h->Draw("PE");
        else h->Draw("PE same");
      } // iBg
      AddOneLineText(sliceVarLow + " < " + sliceVarName + " < " + sliceVarUp, {0.8, 0.95, 0.9, 0.99});
      leg->Draw("same");
      cc.Print((fileOutName + ".pdf" + priBra).c_str(), "pdf");
      priBra = "";
    } // sliceRanges
    CloseCanvasPrinting({fileOutName});
  };
  PrintCanvases(pTRanges, "pT", 0);
  PrintCanvases(lifetimeRanges, "T", 2);

  const std::string fileOutName = "total_Pgt" + promptBdtCut;
  TCanvas cc("cc", "");
  cc.SetCanvasSize(1200, 800);
  TLegend* leg = new TLegend(0.8, 0.4, 0.95, 0.8);
  leg->SetHeader(("bdt(prompt) > " + promptBdtCut).c_str());
  for (int iBg = 1; iBg <= 10; iBg++) {
    const std::string sBgCut = HelperGeneral::to_string_with_precision((11 - iBg) * 0.01, 2);
    const std::string histoName =
            dataType + "/BGlt" + sBgCut + "/hMass_BGlt" + sBgCut + "_Pgt" + promptBdtCut;
    TH1* h = GetObjectWithNullptrCheck<TH1>(fileIn, histoName);
    h->UseCurrentStyle();
    if (rebinFactor != 1) h->Rebin(rebinFactor);
    h->SetLineColor(palette.at((iBg - 1) % palette.size()));
    leg->AddEntry(h, ("bdt(BG) < " + sBgCut).c_str(), "L");
    if (iBg == 1) h->Draw("PE");
    else h->Draw("PE same");
  } // iBg
  leg->Draw("same");
  cc.Print((fileOutName + ".pdf").c_str(), "pdf");

}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./mass_bdt_qa fileName (dataType=data)" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileName = argv[1];
  const std::string dataType = argc > 2 ? argv[2] : "data";

  if(dataType != "data" && dataType != "prompt" && dataType != "nonprompt") {
    throw std::runtime_error("mass_bdt_qa::main(): dataType must be 'data', 'prompt' or 'nonprompt'");
  }

  MassBdtQa(fileName, dataType);

  return 0;
}