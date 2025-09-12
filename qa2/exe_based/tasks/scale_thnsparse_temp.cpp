//
// Created by oleksii on 12.09.2025.
//

#include "HelperGeneral.hpp"
#include "HelperMath.hpp"

#include <TAxis.h>
#include <THnSparse.h>
#include <TH1.h>

#include <iostream>

using namespace HelperGeneral;

void PerformScaling(THnSparse* histoIn, int nDim, const TH1* histoWeight);

void ScaleTHnSparse(const std::string& fileName, const std::string& fileNamePtWeight) {
  TFile* fileIn = OpenFileWithNullptrCheck(fileName);
  TFile* fileWeight = OpenFileWithNullptrCheck(fileNamePtWeight);
  TFile* fileOut = TFile::Open("scale_thnsparse_temp.root", "recreate");

  THnSparse* histoIn = GetObjectWithNullptrCheck<THnSparse>(fileIn, "hf-task-lc/hnLcVarsGen");
  THnSparse* histoOut = dynamic_cast<THnSparse*>(histoIn->Clone());
  histoOut->SetName((static_cast<std::string>(histoIn->GetName()) + "_W").c_str());
  TH1* histoWeight = GetObjectWithNullptrCheck<TH1>(fileWeight, "histoWeight_pT_0_20");

  PerformScaling(histoOut,0, histoWeight);

  fileOut->cd();
  histoIn->Write();
  histoOut->Write();

  fileOut->Close();
  fileWeight->Close();
  fileIn->Close();
}

void PerformScaling(THnSparse* histoIn, int nDim, const TH1* histoWeight) {
  // Optionally enable Sumw2 (errors), if you want to scale errors properly
  histoIn->Sumw2();

  // Number of dimensions
  int nd = histoIn->GetNdimensions();

  // Buffers for coordinates
  // coords[i] will hold bin index along axis i
  std::vector<int> coords(nd);

  // Loop over filled bins
  Long64_t nfilled = histoIn->GetNbins();  // number of *filled* bins
  for (Long64_t lin = 0; lin < nfilled; ++lin) {
    // Get the current content, and the axis‐coordinates
    double content = histoIn->GetBinContent(lin, coords.data());
    if (content == 0) continue;  // maybe skip zeros, though filled bins usually nonzero

    // Get the value along the axisScale (e.g., Pt axis)
    int binIndex = coords[nDim];  // Note: this is the bin number along that axis
    // Find the physical value: maybe bin center
    TAxis *ax = histoIn->GetAxis(nDim);
    double axisValue = ax->GetBinCenter(binIndex);

    // Compute scale factor
    double sf = InterpolateTH1SuppressWarning(histoWeight, axisValue);

    // Set new content
    double newContent = content * sf;
    histoIn->SetBinContent(coords.data(), newContent);

    // If you want also scale the error (assuming Sumw2 is enabled):
    double err = histoIn->GetBinError(lin);
    double newErr = err * sf;
    histoIn->SetBinError(coords.data(), newErr);
  }
}

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./scale_thnsparse_temp fileName filePtWeight" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileName = argv[1];
  const std::string fileNamePtWeight = argv[2];

  ScaleTHnSparse(fileName, fileNamePtWeight);

  return 0;
}