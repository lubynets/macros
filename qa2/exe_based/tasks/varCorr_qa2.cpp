//
// Created by oleksii on 14.04.25.
//

#include "HelperGeneral.hpp"
#include "HelperPlot.hpp"

#include <TExec.h>
#include <TFile.h>
#include <TH2.h>
#include <TStyle.h>

#include <iostream>
#include <string>

using namespace HelperGeneral;
using namespace HelperPlot;

void VarCorrQa2(const std::string& fileName, const std::string& mcOrData) {
  const std::vector<std::string> vars {"nSigTpcPr", "nSigTpcPi", "nSigTpcKa", "ldl", "chi2Topo", "chi2PrimPr", "chi2PrimPi", "chi2PrimKa", "chi2Geo", "mass"};
  std::vector<std::string> dataTypes;
  if(mcOrData == "mc") {
    dataTypes.emplace_back("prompt");
    dataTypes.emplace_back("nonPrompt");
  } else if (mcOrData == "data") {
    dataTypes.emplace_back("background");
  }
  const std::vector<std::string> ptCutEdges{"0", "2", "5", "8", "12", "20"};
  bool setLogz{true};

  TFile* fileIn = OpenFileWithNullptrCheck(fileName);

  for(const auto& dt : dataTypes) {
    for(int iPt=0, nPts=ptCutEdges.size()-1; iPt<nPts; iPt++) {
      const std::string fileOutName = dt + ".pT_" + ptCutEdges.at(iPt) + "_" + ptCutEdges.at(iPt+1);
      std::string printingBracket = "(";

      LoadMacro("styles/varCorr.style.cc");
      TH2F* histoCorr = new TH2F("histoCorr", "", vars.size(), 0, vars.size(), vars.size(), 0, vars.size());
      for(int iVar=0, nVars=vars.size(); iVar<nVars; iVar++) {
        histoCorr->GetXaxis()->SetBinLabel(iVar+1, vars.at(nVars-1-iVar).c_str());
        histoCorr->GetYaxis()->SetBinLabel(iVar+1, vars.at(nVars-1-iVar).c_str());
        histoCorr->SetBinContent(iVar+1, iVar+1, 1.f);
      }
      histoCorr->GetZaxis()->SetRangeUser(-1., 1.);
      histoCorr->GetXaxis()->SetLabelOffset(0.01);

      for(int iVar=0, nVars=vars.size(); iVar<nVars; iVar++) {
        const std::string histoName1D = dt + "/pT_" + ptCutEdges.at(iPt) + "_" + ptCutEdges.at(iPt+1) + "/" + vars.at(iVar);
        TH1* histo1D = GetObjectWithNullptrCheck<TH1>(fileIn, histoName1D);
        LoadMacro("styles/mc_qa2.style.cc");
        histo1D->UseCurrentStyle();
        HistoQuantities histo1D_stat = EvaluateHistoQuantities(histo1D);
        TPaveText histo1D_stat_text = ConvertHistoQuantitiesToText(histo1D_stat, 0.70, 0.6, 0.90, 0.8);
        TCanvas c1("c1", "");
        c1.SetCanvasSize(1200, 800);
        histo1D->Draw();
        histo1D_stat_text.Draw("same");
        c1.Print((fileOutName + ".1D.pdf" + printingBracket).c_str(), "pdf");
        for(int jVar=iVar+1; jVar<nVars; jVar++) {
          const std::string histoName2D = dt + "/pT_" + ptCutEdges.at(iPt) + "_" + ptCutEdges.at(iPt+1) + "/" + vars.at(iVar) + "_vs_" + vars.at(jVar);
          const std::string histoNameInv2D = dt + "/pT_" + ptCutEdges.at(iPt) + "_" + ptCutEdges.at(iPt+1) + "/" + vars.at(jVar) + "_vs_" + vars.at(iVar);
          TH2* histo2D;
          try {
            histo2D = GetObjectWithNullptrCheck<TH2>(fileIn, histoName2D);
          } catch (std::exception&) {
            histo2D = GetObjectWithNullptrCheck<TH2>(fileIn, histoNameInv2D);
          }
          LoadMacro("styles/varCorr.style.cc");
          histo2D->UseCurrentStyle();
          TCanvas cc("cc", "cc", 1200, 800);
          cc.SetLogz(setLogz);
          gStyle->SetPalette(kBird);
          histo2D->Draw("colz");
          AddOneLineText(dt, {0.8, 0.95, 0.9, 0.99});
          AddOneLineText("#it{p}_{T}#in (" + ptCutEdges.at(iPt) + "; " + ptCutEdges.at(iPt+1) + ") GeV/#it{c}", {0.2, 0.95, 0.4, 0.99});
          const double corrCoef = histo2D->GetCorrelationFactor();
          const double errCorrCoef = 1./std::sqrt(histo2D->GetEffectiveEntries());
          histoCorr->SetBinContent(nVars-iVar, nVars-jVar, corrCoef);
          histoCorr->SetBinContent(nVars-jVar, nVars-iVar, corrCoef);
          AddOneLineText("r = (" + to_string_with_significant_figures(corrCoef*100, 2) + " #pm " + to_string_with_significant_figures(errCorrCoef*100, 2) + ")%", {0.6, 0.95, 0.7, 0.99});
          cc.Print((fileOutName + ".pdf" + printingBracket).c_str(), "pdf");
          printingBracket = "";
        } // jVar
      } // iVar
      CloseCanvasPrinting({fileOutName + ".1D"});
      TCanvas cc("cc", "cc");
      cc.SetCanvasSize(1200, 1000);
      cc.SetRightMargin(0.14);
      cc.SetLeftMargin(0.16);
      cc.SetBottomMargin(0.10);
      cc.SetTopMargin(1. - cc.GetBottomMargin() - (1.-cc.GetLeftMargin()-cc.GetRightMargin())*cc.GetWw()/cc.GetWh());
      gStyle->SetPalette(kLightTemperature);
      histoCorr->Draw("colz");
      AddOneLineText(dt, {0.8, 0.95, 0.9, 0.99});
      AddOneLineText("#it{p}_{T}#in (" + ptCutEdges.at(iPt) + "; " + ptCutEdges.at(iPt+1) + ") GeV/#it{c}", {0.2, 0.95, 0.4, 0.99});
      cc.Print((fileOutName + ".pdf)").c_str(), "pdf");

      delete histoCorr;
    } // ptCutEdges
  } // dataTypes

}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./varCorr_qa2 fileName mcOrData" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileName = argv[1];
  const std::string mcOrData = argv[2];
  if(mcOrData != "mc" && mcOrData != "data") throw std::runtime_error("varCorr_qa::main(): mcOrData must be either 'mc' or 'data'");

  VarCorrQa2(fileName, mcOrData);

  return 0;
}
