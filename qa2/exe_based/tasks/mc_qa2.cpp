#include "Helper.hpp"

#include <TPaveText.h>
#include <TFile.h>
#include <TROOT.h>
#include <TH2.h>
#include <TCanvas.h>

#include <iostream>

using namespace Helper;

void mc_qa2(const std::string& fileName, int prompt_or_nonprompt) {
  TString currentMacroPath = __FILE__;
  TString directory = currentMacroPath(0, currentMacroPath.Last('/'));
  gROOT->Macro( directory + "/../styles/mc_qa2.style.cc" );

  if(prompt_or_nonprompt !=1 && prompt_or_nonprompt != 2) {
    throw std::runtime_error("prompt_or_nonprompt must be 1 or 2");
  }

  const std::string promptness = prompt_or_nonprompt == 1 ? "prompt" : "nonprompt";

  TFile* fileIn = TFile::Open(fileName.c_str());
  if(fileIn == nullptr) {
    throw std::runtime_error("fileIn == nullptr");
  }

  std::string fileOutName = "mc_qa2";

  struct Variable {
    std::string name_;
    bool log_mc_;
    bool log_rec_;
    bool log_res_;
    bool log_corr_;
    bool log_pull_;
  };

  std::vector<Variable> vars {
//  name    logmc  logrec logres logcorr logpull
    {"P",   false, false, false, true, false},
    {"Pt",  false, false, false, true, false},
    {"Xsv", false, false, false, true, false},
    {"Ysv", false, false, false, true, false},
    {"Zsv", false, false, false, true, false},
    {"L",   false, true,  true,  true, true},
    {"T",   false, true,  true,  true, true},
  };

  std::string printing_bracket;
  for(int iVar=0; iVar<vars.size(); iVar++) {
    if(iVar == 0) printing_bracket = "(";
    else if(iVar==vars.size()-1) printing_bracket = ")";
    else printing_bracket = "";

    auto PrintCanvas1D = [&] (const std::string& histoType, const std::string& dirName, bool logy, const std::string& oneLineText = "") {
      TH1D* h = fileIn->Get<TH1D>(("PullsAndResiduals/" + dirName + "_" + promptness + "_total/" + histoType + "_" + vars.at(iVar).name_).c_str());
      if(h == nullptr) {
        throw std::runtime_error(("h == nullptr for " + histoType + " " + vars.at(iVar).name_).c_str());
      }
      h->UseCurrentStyle();
      TCanvas cc("cc", "cc", 1200, 800);
      cc.SetLogy(logy);
      h->Draw();
      AddOneLineText(promptness, 0.74, 0.86, 0.87, 0.90);
      AddOneLineText(oneLineText, 0.74, 0.78, 0.87, 0.84);
      HistoQuantities quant = EvaluateHistoQuantities(h);
      TPaveText quant_text = ConvertHistoQuantitiesToText(quant, 0.70, 0.6, 0.90, 0.8);
      quant_text.Draw("same");
      cc.Print((fileOutName + "_" + histoType + ".pdf" + printing_bracket).c_str(), "pdf");
    };

    PrintCanvas1D("mc", "Simulated", vars.at(iVar).log_mc_, "MC matched with reco only");
    PrintCanvas1D("rec", "Candidates", vars.at(iVar).log_rec_, "Rec matched to MC-true only");
    PrintCanvas1D("res", "Candidates_Simulated", vars.at(iVar).log_res_);
    PrintCanvas1D("pull", "Candidates_Simulated", vars.at(iVar).log_pull_);

    TH2D* hcorr = fileIn->Get<TH2D>(("PullsAndResiduals/Candidates_Simulated_"  + promptness + "_total/corr_" + vars.at(iVar).name_).c_str());
    if(hcorr == nullptr) {
      throw std::runtime_error(("hcorr == nullptr for " + vars.at(iVar).name_).c_str());
    }
    hcorr->UseCurrentStyle();

    TCanvas ccCorr("ccCorr", "ccCorr", 1200, 800);
    ccCorr.SetRightMargin(0.16);
    ccCorr.SetLogz(vars.at(iVar).log_corr_);
    hcorr->Draw("colz");
    AddOneLineText(promptness, 0.74, 0.86, 0.87, 0.90);

    ccCorr.Print((fileOutName + "_corr.pdf" + printing_bracket).c_str(), "pdf");
  } // vars
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./mc_qa2 fileName (prompt_or_nonprompt)" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileName = argv[1];
  const int prompt_or_nonprompt = argc>2 ? atoi(argv[2]) : 1;

  mc_qa2(fileName, prompt_or_nonprompt);

  return 0;
}
