//
// Created by oleksii on 29.04.25.
//
#include "BdtEfficiencyCalculator.hpp"

#include "HelperGeneral.hpp"
#include "HelperPlot.hpp"

#include <TFile.h>
#include <TH2.h>

#include <iostream>

using namespace HelperGeneral;
using namespace HelperPlot;

void BdtEfficiency(const std::string& fileNameMc, const std::string& fileNameData, const std::string& xCutDirection, const std::string& yCutDirection) {
  LoadMacro("styles/bdt_qa.style.cc");
  TFile* fileInMc = OpenFileWithNullptrCheck(fileNameMc);
  TFile* fileInData = OpenFileWithNullptrCheck(fileNameData);

  struct DataType {
    std::string name_;
    bool is_mc_;
    bool is_oriented_; // whether we want to maximize (true) its efficiency or minimize (false)
  };

  std::vector<DataType> dataTypes {
    {"prompt", true, true},
    {"nonPrompt", true, false},
    {"background", false, false},
  };

  const std::string fileOutYieldName = "bdt_yield";
  const std::string fileOutEffName = "bdt_eff";

  std::vector<BdtEfficiencyCalculator*> becs(dataTypes.size(), nullptr);

  std::string priBra = "(";
  int iDt{0};
  for(const auto& dt : dataTypes) {
    TH2* histoBdt{nullptr};
    if(dt.is_mc_) histoBdt = GetObjectWithNullptrCheck<TH2>(fileInMc, dt.name_ + "/hBdt2D");
    else          histoBdt = GetObjectWithNullptrCheck<TH2>(fileInData, dt.name_ + "/hBdt2D");
    histoBdt->UseCurrentStyle();
    histoBdt->GetZaxis()->SetTitleOffset(1.35);

    becs.at(iDt) = new BdtEfficiencyCalculator(histoBdt, xCutDirection, yCutDirection);
    becs.at(iDt)->Run();
    TH2D* histoBdtOut = dt.is_oriented_ ? becs.at(iDt)->GetEfficiencyHistogram() : becs.at(iDt)->GetRejectionHistogram();
    TCanvas ccEff("ccEff", "");
    TCanvas ccYield("ccYield", "");
    ccYield.SetLogz();
    for(const auto& cc : {&ccEff, &ccYield}) {
      cc->SetCanvasSize(1000, 1000);
    }
    ccYield.cd();
    histoBdt->Draw("colz");
    AddOneLineText(dt.name_, {0.45, 0.95, 0.55, 0.99});
    ccEff.cd();
    histoBdtOut->Draw("colz");
    AddOneLineText(dt.name_, {0.45, 0.95, 0.55, 0.99});
    ccYield.Print((fileOutYieldName + ".pdf" + priBra).c_str(), "pdf");
    ccEff.Print((fileOutEffName + ".pdf" + priBra).c_str(), "pdf");
    priBra = "";
    ++iDt;
  } // dataTypes
  CloseCanvasPrinting({fileOutYieldName, fileOutEffName});

  TH2D* histoSignalEffRatio = dynamic_cast<TH2D*>(becs.at(1)->GetEfficiencyHistogram()->Clone()); // nonPrompt efficiency histo
  histoSignalEffRatio->Divide(becs.at(0)->GetEfficiencyHistogram()); // divide by prompt efficiency
  histoSignalEffRatio->GetZaxis()->SetTitle("#varepsilon_{nonPrompt} / #varepsilon_{prompt}");
  histoSignalEffRatio->GetZaxis()->SetRangeUser(0.01, 1);
  TCanvas ccRatio("ccRatio", "");
  ccRatio.SetLogz();
  ccRatio.SetCanvasSize(1000, 1000);
  histoSignalEffRatio->Draw("colz");
  ccRatio.Print("bdt_signal_ratio.pdf", "pdf");
}


int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./bdt_efficiency fileNameMc fileNameData (xCutDirection=up yCutDirection=low)" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileNameMc = argv[1];
  const std::string fileNameData = argv[2];
  const std::string xCutDirection = argc > 3 ? argv[3] : "up";
  const std::string yCutDirection = argc > 4 ? argv[4] : "low";

  BdtEfficiency(fileNameMc, fileNameData, xCutDirection, yCutDirection);

  return 0;
}