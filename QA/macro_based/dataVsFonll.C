struct DataPoint {
  double cross_section_;
  double stat_error_;
  double syst_error_up_;
  double syst_error_down_;
};

double EvalErrorDataPoint(const DataPoint& datapoint, bool isIncludeSysError=false);

void dataVsFonll() {
  gROOT->Macro("/home/oleksii/alidir/macros_on_git/qa2/exe_based/styles/mc_qa2.style.cc");

  const double picoToMicro = 1.e-6;

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

#include "fonll.h"

  TGraphErrors* grData = new TGraphErrors();
  grData->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  grData->GetYaxis()->SetTitle("d^{2}#sigma / d#it{p}_{T}dy (#mub#upoint#it{c}/GeV)");
  grData->SetLineColor(kBlue);
  grData->SetMarkerColor(kBlue);
  for(int iPoint=0; iPoint<nPoints; ++iPoint) {
    grData->SetPoint(iPoint, (pTEdges.at(iPoint) + pTEdges.at(iPoint+1))/2, sigmas.at(iPoint).cross_section_);
    grData->SetPointError(iPoint, 0., EvalErrorDataPoint(sigmas.at(iPoint)));
  }

  grData->GetYaxis()->SetRangeUser(grData->GetYaxis()->GetXmin(), grData->GetYaxis()->GetXmax()*2);

  TGraphErrors* grFonll = new TGraphErrors();
  grFonll->SetLineColor(kRed);
  grFonll->SetMarkerColor(kRed);
  for(int iPoint=0, nPoints=sigmas_fonll.size(); iPoint<nPoints; ++iPoint) {
    grFonll->SetPoint(iPoint, sigmas_fonll.at(iPoint).first, sigmas_fonll.at(iPoint).second*picoToMicro);
  }

  TLegend* leg = new TLegend(0.75, 0.6, 0.95, 0.8);
  leg->AddEntry(grData, "data", "PE");
  leg->AddEntry(grFonll, "FONLL", "L");

  TCanvas* cc = new TCanvas("cc", "");
  cc->SetCanvasSize(1200, 800);
  cc->SetLogy();

  grData->Draw("APE");
  grFonll->Draw("L same");
  leg->Draw("same");

  cc->Print("dataVsFonll.pdf", "pdf");
}

double EvalErrorDataPoint(const DataPoint& datapoint, bool isIncludeSysError) {
  const double syst_error_ave = isIncludeSysError ? (datapoint.syst_error_up_ + datapoint.syst_error_down_) / 2. : 0.;
  return std::sqrt(datapoint.stat_error_*datapoint.stat_error_ + syst_error_ave*syst_error_ave);
}
