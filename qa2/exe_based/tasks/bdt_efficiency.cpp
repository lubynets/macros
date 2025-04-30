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
    for(int iBin=1; iBin<=x_n_bins_; iBin++) {
      for(int jBin=1; jBin<=y_n_bins_; jBin++) {
        const double eff = GetEfficiency(iBin, jBin);
        histo_eff_->SetBinContent(iBin, jBin, eff);
        histo_rej_->SetBinContent(iBin, jBin, 1. - eff);
      } // iBin
    } // jBin
    for(const auto& h : {histo_eff_, histo_rej_}) {
      h->Scale(100.);
      h->SetMinimum(0);
      h->SetMaximum(100);
    }
  }

  TH2D* GetEfficiencyHistogram() const { return histo_eff_; }

  TH2D* GetRejectionHistogram() const { return histo_rej_; }

 private:
  void Init() {
    x_n_bins_ = histo_in_->GetNbinsX();
    y_n_bins_ = histo_in_->GetNbinsY();
    integral_ = histo_in_->Integral(1, x_n_bins_, 1, y_n_bins_);

    histo_eff_ = dynamic_cast<TH2D*>(histo_in_->Clone());
    histo_eff_->Reset();
    std::string xTitle = histo_in_->GetXaxis()->GetTitle();
    if(x_is_upper_cut_) xTitle += " upper";
    else                xTitle += " lower";
    xTitle += " value";
    histo_eff_->GetXaxis()->SetTitle(xTitle.c_str());
    std::string yTitle = histo_in_->GetYaxis()->GetTitle();
    if(y_is_upper_cut_) yTitle += " upper";
    else                yTitle += " lower";
    yTitle += " value";
    histo_eff_->GetYaxis()->SetTitle(yTitle.c_str());
    histo_eff_->GetZaxis()->SetTitle("#varepsilon, %");

    histo_rej_ = dynamic_cast<TH2D*>(histo_eff_->Clone());
    histo_rej_->GetZaxis()->SetTitle("(1 - #varepsilon), %");
  }

  double GetEfficiency(int binx, int biny) const {
    const int fromX = x_is_upper_cut_ ? 1 : binx+1;
    const int toX = x_is_upper_cut_ ? binx-1 : x_n_bins_;
    const int fromY = y_is_upper_cut_ ? 1 : biny+1;
    const int toY = y_is_upper_cut_ ? biny-1 : y_n_bins_;
    const double integralInternal = (toX > fromX && toY > fromY) ? histo_in_->Integral(fromX, toX, fromY, toY) : 0.;
    const double integralEdgeVert = toY > fromY ? histo_in_->Integral(binx, binx, fromY, toY) / 2 : 0.;
    const double integralEdgeHoriz = toX > fromX ? histo_in_->Integral(fromX, toX, biny, biny) / 2 : 0.;
    const double integralCorner = histo_in_->GetBinContent(binx, biny) / 4;
    return (integralInternal + integralEdgeVert + integralEdgeHoriz + integralCorner) / integral_;
  }

  const TH2* histo_in_;
  TH2D* histo_eff_;
  TH2D* histo_rej_;
  bool x_is_upper_cut_;
  bool y_is_upper_cut_;
  int x_n_bins_;
  int y_n_bins_;
  double integral_;
};

void BdtEfficiency(const std::string& fileNameMc, const std::string& fileNameData, const std::string& xCutDirection, const std::string& yCutDirection) {
  LoadMacro("styles/bdt_qa.style.cc");
  TFile* fileInMc = OpenFileWithNullptrCheck(fileNameMc);
  TFile* fileInData = OpenFileWithNullptrCheck(fileNameData);

  struct DataType {
    std::string name_;
    bool is_mc_;
    bool is_oriented_; // whether we want to maximize its efficiency or minimize
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
    if(dt.is_mc_) histoBdt = GetObjectWithNullptrCheck<TH2>(fileInMc, dt.name_ + "/hBdt");
    else          histoBdt = GetObjectWithNullptrCheck<TH2>(fileInData, dt.name_ + "/hBdt");
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