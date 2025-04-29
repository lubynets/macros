//
// Created by oleksii on 29.04.25.
//
#include "Helper.hpp"

#include <TFile.h>
#include <TH2.h>

#include <iostream>

using namespace Helper;

class BdtEfficiencyCalculator {
 public:
  BdtEfficiencyCalculator() = delete;
  BdtEfficiencyCalculator(const TH2* histo, const std::string& xCutDirection, const std::string& yCutDirection): histo_in_(histo) {
    if(histo == nullptr) throw std::runtime_error("BdtEfficiencyCalculator::BdtEfficiencyCalculator(): histo == nullptr");

    if(xCutDirection == "up") x_is_upper_cut_ = true;
    else if(xCutDirection == "low") x_is_upper_cut_ = false;
    else throw std::runtime_error("BdtEfficiencyCalculator::BdtEfficiencyCalculator(): xCutDirection must be either 'up' or 'low'");

    if(yCutDirection == "up") y_is_upper_cut_ = true;
    else if(yCutDirection == "low") y_is_upper_cut_ = false;
    else throw std::runtime_error("BdtEfficiencyCalculator::BdtEfficiencyCalculator(): yCutDirection must be either 'up' or 'low'");

    Init();
  }

  void Run() {
    for(int iBin=1; iBin<x_n_bins_; iBin++) {
      for(int jBin=1; jBin<y_n_bins_; jBin++) {
        histo_out_->SetBinContent(iBin, jBin, GetEfficiency(iBin, jBin));
      } // iBin
    } // jBin
  }

  TH2D* GetEfficiencyHistogram() const { return histo_out_; }

 private:
  void Init() {
    x_n_bins_ = histo_in_->GetNbinsX();
    y_n_bins_ = histo_in_->GetNbinsY();
    integral_ = histo_in_->Integral(1, x_n_bins_, 1, y_n_bins_);

    histo_out_ = dynamic_cast<TH2D*>(histo_in_->Clone());
    histo_out_->Reset();
    std::string xTitle = histo_in_->GetXaxis()->GetTitle();
    if(x_is_upper_cut_) xTitle += " upper";
    else                xTitle += " lower";
    xTitle += " value";
    histo_out_->GetXaxis()->SetTitle(xTitle.c_str());
    std::string yTitle = histo_in_->GetYaxis()->GetTitle();
    if(y_is_upper_cut_) yTitle += " upper";
    else                yTitle += " lower";
    yTitle += " value";
    histo_out_->GetYaxis()->SetTitle(yTitle.c_str());
  }

  double GetEfficiency(int binx, int biny) const {
    const int fromX = x_is_upper_cut_ ? 1 : binx+1;
    const int toX = x_is_upper_cut_ ? binx-1 : x_n_bins_;
    const int fromY = y_is_upper_cut_ ? 1 : biny+1;
    const int toY = y_is_upper_cut_ ? biny-1 : y_n_bins_;
    const double integralInternal = histo_in_->Integral(fromX, toX, fromY, toY);
    const double integralEdge = histo_in_->Integral(fromX, toX, biny, biny) / 2 +
                                histo_in_->Integral(binx, binx, fromY, toY) / 2 +
                                histo_in_->GetBinContent(binx, biny) / 4;
    return (integralInternal + integralEdge) / integral_;
  }

  const TH2* histo_in_;
  TH2D* histo_out_;
  bool x_is_upper_cut_;
  bool y_is_upper_cut_;
  int x_n_bins_;
  int y_n_bins_;
  double integral_;
};

void BdtEfficiency(const std::string& fileNameMc, const std::string& fileNameData, const std::string& xCutDirection, const std::string& yCutDirection) {
  LoadMacro("styles/mc_qa2.style.cc");
  TFile* fileInMc = OpenFileWithNullptrCheck(fileNameMc);
  TFile* fileInData = OpenFileWithNullptrCheck(fileNameData);

  std::vector<std::string> dataTypes{"prompt", "nonPrompt", "background"};

  const std::string fileOutName = "bdt_efficiency.out";

  std::string priBra = "(";
  for(const auto& dt : dataTypes) {
    TH2* histoBdt{nullptr};
    if(dt == "prompt" || dt == "nonPrompt") histoBdt = GetObjectWithNullptrCheck<TH2>(fileInMc, dt + "/hBdt");
    else                                    histoBdt = GetObjectWithNullptrCheck<TH2>(fileInData, dt + "/hBdt");

    BdtEfficiencyCalculator bec(histoBdt, xCutDirection, yCutDirection);
    bec.Run();
    TH2D* histoBdtEff = bec.GetEfficiencyHistogram();
    TCanvas cc("cc", "");
    cc.SetCanvasSize(1000, 1000);
    histoBdtEff->Draw("colz");
    AddOneLineText(dt, {0.8, 0.95, 0.9, 0.99});
    cc.Print((fileOutName + ".pdf" + priBra).c_str(), "pdf");
    priBra = "";
  }
  CloseCanvasPrinting({fileOutName});
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