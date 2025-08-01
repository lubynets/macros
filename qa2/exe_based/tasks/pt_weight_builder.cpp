//
// Created by oleksii on 30.07.25.
//
#include "HelperGeneral.hpp"
#include "HelperMath.hpp"

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TLegend.h>

#include <string>

using namespace HelperGeneral;
using namespace HelperMath;

void pt_weight_builder(const std::string& fileNamePtGen, const std::string& fileNamePtFit, bool isGenHistoAccumulated) {
  LoadMacro("styles/mc_qa2.style.cc");

  const std::string fileOutName = "ptWeight";

  TFile* fileGen = OpenFileWithNullptrCheck(fileNamePtGen, "read");
  TFile* fileFit = OpenFileWithNullptrCheck(fileNamePtFit, "read");

  TH1* histoGen = GetObjectWithNullptrCheck<TH1>(fileGen, "histoPtGen");
  TF1* funcFit = GetObjectWithNullptrCheck<TF1>(fileFit, "tsallisFit");
  histoGen->UseCurrentStyle();
  funcFit->UseCurrentStyle();

  TFile* fileOut = TFile::Open((fileOutName + ".root").c_str(), "recreate");
  histoGen->Write();

  auto ProcessPtRange = [&] (double loPt, double hiPt) {
    const std::string ptRangeString = "pT_" + HelperGeneral::to_string_with_precision(loPt, 0) + "_" + HelperGeneral::to_string_with_precision(hiPt, 0);

    TH1* histoCut = CutSubHistogram(histoGen, loPt, hiPt);

    if(isGenHistoAccumulated) histoCut->Scale(1., "width");

    if(loPt < funcFit->GetXmin() || hiPt > funcFit->GetXmax()) throw std::runtime_error("loPt < funcFit->GetXmin() || hiPt > funcFit->GetXmax()");
  
    const double integralGen = histoCut->Integral("width");
    const double integralFit = funcFit->Integral(loPt, hiPt);
  
    TH1* histoCutNorm = dynamic_cast<TH1*>(histoCut->Clone());
    histoCutNorm->GetYaxis()->SetTitle("d^{2}#sigma / d#it{p}_{T}dy (a.u.)");
    Sumw2IfNotYet(histoCutNorm);
    histoCutNorm->SetMarkerStyle(0);
    histoCutNorm->SetLineWidth(2);
  
    TF1* funcFitNorm = dynamic_cast<TF1*>(funcFit->Clone());
    histoCutNorm->Scale(1./integralGen);
    funcFitNorm->SetParameter(0, funcFitNorm->GetParameter(0) / integralFit); // FIXME function implementation-defined, relies on the [0] parameter to be a common factor

    TH1* histoWeight = dynamic_cast<TH1*>( histoCutNorm->Clone());
    histoWeight->GetYaxis()->SetTitle("weight, Tsallis / Pythia");
    histoWeight->GetYaxis()->SetRangeUser(0.5, 2);
    DivideHistoByFunction(histoWeight, funcFitNorm, "I");
    InvertHisto(histoWeight);
  
    TCanvas ccShapes("ccShapes", "");
    ccShapes.SetCanvasSize(1200, 800);
    ccShapes.SetLogy();
    histoCutNorm->Draw("HIST");
    funcFitNorm->Draw("same");
    TLegend leg(0.7, 0.82, 0.9, 0.9);
    leg.AddEntry(histoCutNorm, "Pythia", "L");
    leg.AddEntry(funcFitNorm, "Tsallis", "L");
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
  
    histoCutNorm->Write(("histoPtGenNorm_" + ptRangeString).c_str());
    funcFitNorm->Write(("tsallisFitNorm_" + ptRangeString).c_str());
    histoWeight->Write(("histoWeight_" + ptRangeString).c_str());
  };

  ProcessPtRange(0., 20.);

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