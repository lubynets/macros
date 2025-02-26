#include "Helper.hpp"

#include <TPaveText.h>
#include <TFile.h>
#include <TROOT.h>
#include <TH2.h>
#include <TCanvas.h>

#include <iostream>

using namespace Helper;

void mc_qa2(const std::string& fileName) {
  TString currentMacroPath = __FILE__;
  TString directory = currentMacroPath(0, currentMacroPath.Last('/'));
  gROOT->Macro( directory + "/../styles/mc_qa2.style.cc" );

  TFile* fileIn = TFile::Open(fileName.c_str());
  if(fileIn == nullptr) throw std::runtime_error("fileIn == nullptr");

  const std::string fileOutNamePrefix = "mc_qa2";

  struct Variable {
    std::string name_;
    bool log_mc_;
    bool log_rec_;
    bool log_res_;
    bool log_corr_;
    bool log_pull_;
  };

  auto PrintCanvas1D = [&] (const Variable& var,
                            const std::string& histoType,
                            const std::string& dirName,
                            bool logy,
                            const std::string& printing_bracket,
                            const std::string fileOutName,
                            const std::vector<std::string> oneLineTexts={}) {
    TH1D* h = fileIn->Get<TH1D>((dirName + histoType + "_" + var.name_).c_str());
    if(h == nullptr) throw std::runtime_error(("h == nullptr for " + histoType + " " + var.name_).c_str());
    h->UseCurrentStyle();
    TCanvas cc("cc", "cc", 1200, 800);
    cc.SetLogy(logy);
    h->Draw();
    for(int iOLT=0; iOLT<oneLineTexts.size(); iOLT++) {
      AddOneLineText(oneLineTexts.at(iOLT), 0.74, 0.86 - iOLT*0.04, 0.87, 0.90 - iOLT*0.04);
    }
    HistoQuantities quant = EvaluateHistoQuantities(h);
    TPaveText quant_text = ConvertHistoQuantitiesToText(quant, 0.70, 0.6, 0.90, 0.8);
    quant_text.Draw("same");
    cc.Print((fileOutName + "_" + histoType + ".pdf" + printing_bracket).c_str(), "pdf");
  };

  auto PrintCanvasCorr = [&] (const Variable& var,
                              const std::string& dirName,
                              const std::string& printing_bracket,
                              const std::string fileOutName,
                              const std::vector<std::string> oneLineTexts={}) {
    TH2D* hcorr = fileIn->Get<TH2D>((dirName + "corr_" + var.name_).c_str());
    if(hcorr == nullptr) throw std::runtime_error(("hcorr == nullptr for " + var.name_).c_str());
    hcorr->UseCurrentStyle();
    TCanvas ccCorr("ccCorr", "ccCorr", 1200, 800);
    ccCorr.SetRightMargin(0.16);
    ccCorr.SetLogz(var.log_corr_);
    hcorr->Draw("colz");
    for(int iOLT=0; iOLT<oneLineTexts.size(); iOLT++) {
      AddOneLineText(oneLineTexts.at(iOLT), 0.64, 0.86 - iOLT*0.04, 0.77, 0.90 - iOLT*0.04);
    }
    ccCorr.Print((fileOutName + "_corr.pdf" + printing_bracket).c_str(), "pdf");
  };

  std::vector<Variable> varsCand {
    {"P",   false, false, false, true, false},
    {"Pt",  false, false, false, true, false},
    {"Xsv", false, false, false, true, false},
    {"Ysv", false, false, false, true, false},
    {"Zsv", false, false, false, true, false},
    {"L",   false, true,  true,  true, true},
    {"T",   false, true,  true,  true, true},
  };

  std::vector<Variable> varsEve {
    {"Xpv", false, false, false, true, false},
    {"Ypv", false, false, false, true, false},
    {"Zpv", false, false, false, true, false},
  };

  for(auto& promptness : std::vector<std::string>{"prompt", "nonprompt"}) {
    for(int iVar=0; iVar < varsCand.size(); iVar++) {
      const std::string printing_bracket = iVar == 0 ? "(" : iVar == varsCand.size() - 1 ? ")" : "";

      const std::string dirNamePrefix = "PullsAndResiduals/CandidatesQA/";
      const std::string dirNameSuffix = promptness + "_total/";
      const std::string fileOutName = fileOutNamePrefix + "_" + promptness;

      PrintCanvas1D(varsCand.at(iVar), "mc", dirNamePrefix + "Simulated_" + dirNameSuffix, varsCand.at(iVar).log_mc_, printing_bracket, fileOutName, {promptness, "MC matched with reco only"});
      PrintCanvas1D(varsCand.at(iVar), "rec", dirNamePrefix + "Candidates_" + dirNameSuffix, varsCand.at(iVar).log_rec_, printing_bracket, fileOutName, {promptness, "Rec matched to MC-true only"});
      PrintCanvas1D(varsCand.at(iVar), "res", dirNamePrefix + "Candidates_Simulated_" + dirNameSuffix, varsCand.at(iVar).log_res_, printing_bracket, fileOutName, {promptness});
      PrintCanvas1D(varsCand.at(iVar), "pull", dirNamePrefix + "Candidates_Simulated_" + dirNameSuffix, varsCand.at(iVar).log_pull_, printing_bracket, fileOutName, {promptness});

      PrintCanvasCorr(varsCand.at(iVar), dirNamePrefix + "Candidates_Simulated_" + dirNameSuffix, printing_bracket, fileOutName, {promptness});
    } // varsCand
  } // promptnesses

  for(int iVar=0; iVar<varsEve.size(); iVar++) {
    const std::string printing_bracket = iVar == 0 ? "(" : iVar == varsEve.size() - 1 ? ")" : "";

    const std::string dirName = "PullsAndResiduals/EventsQA/Events/";
    const std::string fileOutName = fileOutNamePrefix + "_events";

    PrintCanvas1D(varsEve.at(iVar), "mc", dirName, varsEve.at(iVar).log_mc_, printing_bracket, fileOutName);
    PrintCanvas1D(varsEve.at(iVar), "rec", dirName, varsEve.at(iVar).log_rec_, printing_bracket, fileOutName);
    PrintCanvas1D(varsEve.at(iVar), "res", dirName, varsEve.at(iVar).log_res_, printing_bracket, fileOutName);
    PrintCanvas1D(varsEve.at(iVar), "pull", dirName, varsEve.at(iVar).log_pull_, printing_bracket, fileOutName);

    PrintCanvasCorr(varsEve.at(iVar), dirName, printing_bracket, fileOutName);
  } // varsEve
}

int main(int argc, char* argv[]) {
  if (argc < 1) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./mc_qa2 fileName" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileName = argv[1];

  mc_qa2(fileName);

  return 0;
}
