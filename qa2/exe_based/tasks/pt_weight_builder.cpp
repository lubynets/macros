//
// Created by oleksii on 30.07.25.
//
#include "HelperGeneral.hpp"
#include "HelperMath.hpp"

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <THnSparse.h>
#include <TLegend.h>

#include <string>
#include <string_view>

using namespace HelperGeneral;
using namespace HelperMath;
using namespace std::string_literals;

const std::string_view signalTypeAxisTitle = "candidates type";
const std::string_view pTAxisTitle = "#it{p}_{T}(#Lambda_{c}^{+}) (GeV/#it{c})";

void pt_weight_builder(const std::string& fileNamePtGen, const std::string& fileNamePtFit, bool isGenHistoAccumulated) {
  LoadMacro("styles/mc_qa2.style.cc");

  const double loPt{0.};
  const double hiPt{20.};

  const std::string fileOutName = "ptWeight";

  TFile* fileGen = OpenFileWithNullptrCheck(fileNamePtGen, "read");
  TFile* fileFit = OpenFileWithNullptrCheck(fileNamePtFit, "read");

  THnSparse* histoGenTHn = GetObjectWithNullptrCheck<THnSparse>(fileGen, "hf-task-lc/hnLcVarsGen");
  const std::map<std::string_view, int> axesIndices = MapTHnSparseAxesIndices(histoGenTHn);
  SetTHnSparseAxisRanges(histoGenTHn, axesIndices.at(signalTypeAxisTitle), 1., 2.);

  TFile* fileOut = TFile::Open((fileOutName + ".root").c_str(), "recreate");

  TH1* histoGen = histoGenTHn->Projection(axesIndices.at(pTAxisTitle));
  histoGen->UseCurrentStyle();
  histoGen->Write();

  auto ProcessTsallisFunc = [&] (TH1* histoFit, const bool isMean) {
    const std::string ptRangeString = "pT_" + HelperGeneral::to_string_with_precision(loPt, 0) + "_" + HelperGeneral::to_string_with_precision(hiPt, 0);

    histoGen = CutSubHistogram(histoGen, loPt, hiPt);
    histoFit = CutSubHistogram(histoFit, loPt, hiPt);

    if(isGenHistoAccumulated) histoGen->Scale(1., "width");

    const double integralGen = histoGen->Integral("width");
    const double integralFit = histoFit->Integral("width");
  
    TH1* histoGenNorm = dynamic_cast<TH1*>(histoGen->Clone());
    histoGenNorm->GetYaxis()->SetTitle("d^{2}#sigma / d#it{p}_{T}dy (a.u.)");
    Sumw2IfNotYet(histoGenNorm);
    histoGenNorm->SetMarkerStyle(0);
    histoGenNorm->SetLineWidth(2);
  
    TH1* histoFitNorm = dynamic_cast<TH1*>(histoFit->Clone());
    histoGenNorm->Scale(1./integralGen);
    histoFitNorm->Scale(1./integralFit);
    for (int iBin=1, nBins=histoFitNorm->GetNbinsX(); iBin <= nBins; ++iBin) {
      histoFitNorm->SetBinError(iBin, 0.);
    }

    TH1* histoWeight = dynamic_cast<TH1*>( histoGenNorm->Clone());
    histoWeight->SetName("histoWeight");
    histoWeight->GetYaxis()->SetTitle("weight, Tsallis / Pythia");
    histoWeight->GetYaxis()->SetRangeUser(0.5, 2);
    histoWeight->Sumw2();
    histoWeight->Divide(histoFitNorm);
    InvertHisto(histoWeight);
  
    if(isMean) {
      TCanvas ccShapes("ccShapes", "");
      ccShapes.SetCanvasSize(1200, 800);
      ccShapes.SetLogy();
      histoGenNorm->Draw("HIST");
      histoFitNorm->Draw("same");
      TLegend leg(0.7, 0.82, 0.9, 0.9);
      leg.AddEntry(histoGenNorm, "Pythia", "L");
      leg.AddEntry(histoFitNorm, "Tsallis", "L");
      leg.Draw("same");
      ccShapes.Print((fileOutName + "_" + ptRangeString + ".pdf(").c_str(), "pdf");

      TCanvas ccWeight("ccWeight", "");
      ccWeight.SetCanvasSize(1200, 800);
      histoWeight->Draw();
      TF1 oneline("oneline", "[0]", loPt, hiPt);
      oneline.SetParameter(0, 1);
      oneline.SetLineColor(kBlack);
      oneline.SetLineStyle(7);
      oneline.Draw("same");
      ccWeight.Print((fileOutName + "_" + ptRangeString + ".pdf)").c_str(), "pdf");

      histoGenNorm->Write("histoPtGenNorm");
    }
  
    histoFitNorm->Write(static_cast<TString>(histoFit->GetName()).ReplaceAll("tsallisFitHisto", "tsallisFitHistoNorm"));
    histoWeight->Write(static_cast<TString>(histoFit->GetName()).ReplaceAll("tsallisFitHisto", "histoWeight"));
  };

  TH1* histoFitCent = GetObjectWithNullptrCheck<TH1>(fileFit, "tsallisFitHistocent");
  histoFitCent->UseCurrentStyle();
  ProcessTsallisFunc(histoFitCent, true);
  TH1* histoFitMin = GetObjectWithNullptrCheck<TH1>(fileFit, "tsallisFitHistomin");
  histoFitMin->UseCurrentStyle();
  ProcessTsallisFunc(histoFitMin, false);
  TH1* histoFitMax = GetObjectWithNullptrCheck<TH1>(fileFit, "tsallisFitHistomax");
  histoFitMax->UseCurrentStyle();
  ProcessTsallisFunc(histoFitMax, false);

  fileGen->Close();
  fileFit->Close();
  fileOut->Close();
}

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./pt_weight_builder fileNameGen fileNameFit (isGenHistoAccumulated=true)" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileNameGen = argv[1];
  const std::string fileNameFit = argv[2];
  const bool isGenHistoAccumulated = argc > 3 ? string_to_bool(argv[3]) : true;

  pt_weight_builder(fileNameGen, fileNameFit, isGenHistoAccumulated);

  return 0;
}
