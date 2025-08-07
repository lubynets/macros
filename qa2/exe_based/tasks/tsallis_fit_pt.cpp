//
// Created by oleksii on 28.07.25.
//

#include "HelperGeneral.hpp"
#include "HelperMath.hpp"
#include "HelperPlot.hpp"

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TH1.h>
#include <TMatrixD.h>

#include <array>

using namespace HelperGeneral;
using namespace HelperMath;
using namespace HelperPlot;

struct DataPoint {
    double cross_section_;
    double stat_error_;
    double syst_error_up_;
    double syst_error_down_;
};

double EvalErrorDataPoint(const DataPoint& datapoint, bool isIncludeSysError= true);

void pt_fit(bool isIncludeSysErr, bool isCalculateFitFuncError) {
  LoadMacro("styles/mc_qa2.style.cc");

  // =========================================================================================
  // https://www.hepdata.net/record/ins2697877
  const std::string hepDataName = "www.hepdata.net/record/ins2697877";
  constexpr int nPoints{11};
  const std::array<double, nPoints+1> pTEdges {0., 1., 2., 3., 4., 5., 6., 7., 8., 10., 12., 24.};

  const std::array<DataPoint, nPoints> sigmas {{
    {74.13,  10.237, 8.8693, 8.9308},
    {117.04, 10.146, 11.68,  11.794},
    {72.16,  4.4073, 7.0247, 7.1269},
    {34.96,  1.8675, 3.386,  3.4385},
    {14.834, 0.8375, 1.3392, 1.3781},
    {7.8826, 0.4345, 0.7167, 0.7403},
    {3.4434, 0.2416, 0.3458, 0.3634},
    {1.962,  0.1492, 0.2112, 0.2197},
    {0.8147, 0.0608, 0.0906, 0.0956},
    {0.3155, 0.0395, 0.0368, 0.0388},
    {0.036,  0.0057, 0.0047, 0.0051}
  }};
  // =========================================================================================


  const std::string ptString = "#it{p}_{T}";
  TH1D* histoCrossSec = new TH1D("histoCrossSec", "", nPoints, pTEdges.data());
  for(int iBin=0; iBin<nPoints; ++iBin) {
    histoCrossSec->SetBinContent(iBin+1, sigmas.at(iBin).cross_section_);
    histoCrossSec->SetBinError(iBin+1, EvalErrorDataPoint(sigmas.at(iBin), isIncludeSysErr));
  }
  histoCrossSec->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  histoCrossSec->GetYaxis()->SetTitle("d^{2}#sigma / d#it{p}_{T}dy (#mub#upoint#it{c}/GeV)");
  histoCrossSec->SetLineColor(kBlue);
  histoCrossSec->SetMarkerColor(kBlue);

  TF1* tsallisFit = new TF1("tsallisFit","[0]*x*TMath::Sqrt([1]*[1]+x*x)*TMath::Power(1 + ([2]-1)*TMath::Sqrt([1]*[1]+x*x)/[3],-[2]/([2]-1))", 0., 24.);
  tsallisFit->SetNpx(2400);
  tsallisFit->SetParameters(1.2e+05, 2.28e+00, 1.1e+00, 0.07);
  const std::array<std::string, 4> parNames{"A", "m", "q", "T"};
  for(int iPar=0; iPar<4; ++iPar) {
    tsallisFit->SetParName(iPar, parNames.at(iPar).c_str());
  }
  tsallisFit->SetLineColor(kRed);
  TFitResultPtr tsallisFitResPtr = histoCrossSec->Fit(tsallisFit, "RIS0");
  TMatrixDSym tsallisFitCov = tsallisFitResPtr->GetCovarianceMatrix();

  const float textX1 = 0.62;
  const float textX2 = 0.9;
  const float textY2 = 0.78;
  const float textYStep = 0.04;

  TCanvas ccFit("ccFit", "");
  ccFit.SetCanvasSize(1200, 800);
  ccFit.SetLogy();
  histoCrossSec->Draw();
  tsallisFit->Draw("same");
  auto hepDataText = AddOneLineText("#it{" + hepDataName + "}", {0.3, 0.85, 0.55, 0.90});
  hepDataText->SetTextColor(kGray+2);
  AddOneLineText(parNames.at(0) + ptString + "#sqrt{" + parNames.at(1) + "^{2} + " + ptString + "^{2}} (1 + (" + parNames.at(2) + " - 1)#sqrt{" + parNames.at(1) + "^{2} + " + ptString + "^{2}} / " + parNames.at(3) + ")^{#minus #frac{" + parNames.at(2) + "}{" + parNames.at(2) + " - 1}}", {textX1, textY2 + textYStep, textX2, textY2 + 2*textYStep});
  for(int iPar=0, nPars=tsallisFit->GetNpar(); iPar<nPars; ++iPar) {
    AddOneLineText(static_cast<std::string>(parNames.at(iPar)) + " = " +
                    HelperGeneral::to_string_with_significant_figures(tsallisFit->GetParameter(iPar), 3) + " #pm " +
                    HelperGeneral::to_string_with_significant_figures(tsallisFit->GetParError(iPar), 3), {textX1, textY2 - (iPar+1)*textYStep, textX2, textY2 - iPar*textYStep});
  }
  AddOneLineText("#chi^{2}/ndf = " + HelperGeneral::to_string_with_significant_figures(tsallisFit->GetChisquare()/tsallisFit->GetNDF(), 3), {textX1, textY2 - 5*textYStep, textX2, textY2 - 4*textYStep});
  ccFit.Print("tsallisFitPt.pdf(", "pdf");

  TH1* hRatio = dynamic_cast<TH1*>(histoCrossSec->Clone());
  ScalePlotVertically(hRatio, histoCrossSec, 2);
  hRatio->GetYaxis()->SetTitle("Data / Fit");
  DivideHistoByFunction(hRatio, tsallisFit, "I");

  TCanvas ccRatio("ccRatio", "");
  ScaleCanvasVertically(&ccRatio, &ccFit, 2);
  TF1 oneline("oneline", "[0]", 0, 24);
  oneline.SetParameter(0, 1);
  oneline.SetLineColor(kBlack);
  oneline.SetLineStyle(7);
  hRatio->Draw("HIST PE");
  oneline.Draw("same");
  ccRatio.Print("tsallisFitPt.pdf)", "pdf");

  TFile* fileOut = TFile::Open("pTFit.root", "recreate");
  histoCrossSec->Write("hepData");
  tsallisFit->Write("tsallisFit");
  tsallisFitCov.Write("tsallisFitCov");
  ccFit.Write("ccFit");
  ccRatio.Write("ccRatio");
  fileOut->Close();
}

int main(int argc, char* argv[]) {
  if(argc > 1 && std::strcmp(argv[1], "--help") == 0) {
    std::cout << "./tsallis_pt_fit (isIncludeSysErr=true isCalculateFitFuncError=false)\n";
    return 0;
  }

  const bool isIncludeSysErr = argc > 1 ? HelperGeneral::string_to_bool(argv[1]) : true;
  const bool isCalculateFitFuncError = argc > 2 ? HelperGeneral::string_to_bool(argv[2]) : false;

  pt_fit(isIncludeSysErr, isCalculateFitFuncError);

  return 0;
}

double EvalErrorDataPoint(const DataPoint& datapoint, bool isIncludeSysError) {
  const double syst_error_ave = isIncludeSysError ? (datapoint.syst_error_up_ + datapoint.syst_error_down_) / 2. : 0.;
  return std::sqrt(datapoint.stat_error_*datapoint.stat_error_ + syst_error_ave*syst_error_ave);
}

