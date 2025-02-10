//
// Created by oleksii on 10.02.25.
//
#include "Helper.hpp"

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TString.h>

#include <string>
#include <iostream>
#include <TH2D.h>

using namespace Helper;

void pid_qa2(const std::string& fileName) {

//    const std::string statusDcaFSel = "isDcaFSel";
  const std::string statusDcaFSel = "noDcaFSel";

  const bool setLogY = false;
//  const bool setLogY = true;

  TString currentMacroPath = __FILE__;
  TString directory = currentMacroPath(0, currentMacroPath.Last('/'));
  gROOT->Macro( directory + "/treeKF_qa2.style.cc" );

  TFile* fileIn = TFile::Open(fileName.c_str());
  if(fileIn == nullptr) {
    throw std::runtime_error("fileIn == nullptr");
  }

  const std::string fileOutName = "pid_qa2";

  std::vector<std::string> SignalSpecies {
//    "prompt",
//    "nonprompt",
//    "background",
//    "wrongswap",
    "data",
  };

  std::vector<std::string> PidDetectors{"Tpc", "Tof", "TpcTof"};
  std::vector<std::string> ProngSpecies{"Pr", "Ka", "Pi"};

  for(auto& ss : SignalSpecies) {
    for(auto& prong : ProngSpecies) {
      bool is_first_canvas{true};
      const std::string currentFileOutName = fileOutName + "." + ss + "." + prong + ".pdf";
      for(auto& det : PidDetectors) {
        const std::string histoName = "Candidates_" + ss + "_" + statusDcaFSel + "_total/nSig" + det + prong + "_" + ss + "_" + statusDcaFSel;
        TH1D* histo = fileIn->Get<TH1D>(histoName.c_str());
        if(histo == nullptr) throw std::runtime_error(histoName + " is not present");

        histo->UseCurrentStyle();
        TCanvas cc("cc", "cc", 1200, 800);
        cc.SetLogy(setLogY);
        histo->Draw();
        AddOneLineText(ss, 0.74, 0.82, 0.87, 0.90);
        HistoQuantities quant = EvaluateHistoQuantities(histo);
        TPaveText quant_text = ConvertHistoQuantitiesToText(quant, 0.70, 0.6, 0.90, 0.8);
        quant_text.Draw("same");

        const std::string printing_bracket = is_first_canvas ? "(" : "";
        cc.Print((currentFileOutName + printing_bracket).c_str(), "pdf");
        is_first_canvas = false;
      } // PidDetectors
      const std::string histoName = "Candidates_" + ss + "_" + statusDcaFSel + "_total/nSig2D_" + prong + "_" + ss + "_" + statusDcaFSel;
      TH2D* histo = fileIn->Get<TH2D>(histoName.c_str());
      if(histo == nullptr) throw std::runtime_error(histoName + " is not present");

      histo->UseCurrentStyle();
      TCanvas cc("cc", "cc", 1200, 800);
      histo->Draw();
      AddOneLineText(ss, 0.74, 0.82, 0.87, 0.90);

      cc.Print((currentFileOutName + ")").c_str(), "pdf");
    } // ProngSpecies
  } // SignalSpecies
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./pid_qa2 fileName" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileName = argv[1];

  pid_qa2(fileName);

  return 0;
}
