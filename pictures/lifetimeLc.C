struct DataPoint{
  std::string publication_{};
  std::string collaboration_{};
  int year_{};
  double value_{};
  double stat_err_{};
  double syst_err_{};
};

const std::vector<DataPoint> dataPoints {
  {"",               "ALICE",    2026, -999.,  0.0,  0.0 },
  {"PRL 130 071802", "BELLE II", 2023, 203.20, 0.89, 0.77},
  {"PRD 100 032001", "LHCb",     2019, 202.1,  1.7,  0.9 },
  {"PRL 88 161801",  "FOCUS",    2002, 204.6,  3.4,  2.5 },
  {"PRL 86 5243",    "SELEX",    2001, 198.1,  7.0,  5.6 },
  {"PRL 86 2232",    "CLEO",     2001, 179.6,  6.9,  4.4 },
  {"PRL 70 1755",    "E687",     1993, 215,    16,   8   }
};

const double pdgValue{202.6};
const double pdgError{1.0};

double getCombinedError(const DataPoint& dp) {
  return std::sqrt(dp.stat_err_*dp.stat_err_ + dp.syst_err_*dp.syst_err_);
}

std::string formLabel(const DataPoint& dp) {
  if(!dp.publication_.empty()) {
    return dp.collaboration_ + ", " + dp.publication_ + " (" + std::to_string(dp.year_) + ")";
  } else {
    return dp.collaboration_ + " (" + std::to_string(dp.year_) + ")";
  }
}

void lifetimeLc() {
  const bool isCombineStatSyst{false};

  gStyle->SetCanvasPreferGL(true);
  gStyle->SetPadLeftMargin(0.38);
  gStyle->SetPadRightMargin(0.01);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadTopMargin(0.01);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetFrameLineWidth(4);
  gStyle->SetMarkerSize(2);
  gStyle->SetLineWidth(3);
  gStyle->SetHistLineWidth(3);
  gStyle->SetEndErrorSize(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  gStyle->SetTitleSize(0.06, "X");
  gStyle->SetTitleSize(0.05, "Y");
  gStyle->SetTitleSize(0.05, "Z");

  gStyle->SetLabelSize(0.05, "X");
  gStyle->SetLabelSize(0.06, "Y");
  gStyle->SetLabelSize(0.05, "Z");

  gStyle->SetLabelOffset(0.003, "X");
  gStyle->SetLabelOffset(0.005, "Y");
  gStyle->SetLabelOffset(0.01, "Z");

  gStyle->SetTitleOffset(1.0, "X");
  gStyle->SetTitleOffset(1.46, "Y");
  gStyle->SetTitleOffset(1.15, "Z");

  gStyle->SetNdivisions(208, "xyz");

  gStyle->SetOptStat(0);

  TCanvas* cc = new TCanvas("cc", "", 1200, 800);
  TGraphErrors* grAll = new TGraphErrors();
  TGraphErrors* grSyst = new TGraphErrors();
  grAll->GetXaxis()->SetTitle("#tau_{#Lambda_{c}^{+}} (fs)");
  grAll->SetMarkerColor(kBlack);
  grAll->SetLineColor(kBlack);
  grAll->SetMarkerSize(2);
  grAll->SetMarkerStyle(kFullSquare);
  grSyst->GetXaxis()->SetTitle("#tau_{#Lambda_{c}^{+}} (fs)");
  grSyst->SetMarkerColor(kBlack);
  grSyst->SetLineColor(kBlack);
  grSyst->SetMarkerSize(2);
  grSyst->SetMarkerStyle(kFullSquare);
  grSyst->SetFillStyle(0);

  for(int iDp=0, nDps=dataPoints.size(); iDp<nDps; ++iDp) {
    const auto& dp = dataPoints.at(iDp);
    const auto err1 = isCombineStatSyst ? getCombinedError(dp) : dp.stat_err_;
    const auto errSyst = dp.syst_err_;
    grAll->AddPoint(dp.value_, iDp+0.5);
    grSyst->AddPoint(dp.value_, iDp+0.5);
    grAll->SetPointError(grAll->GetN()-1, err1, 0);
    grSyst->SetPointError(grSyst->GetN()-1, errSyst, 0.15);
  }

  TGraphErrors* grAllThis = dynamic_cast<TGraphErrors*>(grAll->Clone());
  grAllThis->SetMarkerColor(kRed);
  grAllThis->SetLineColor(kRed);
  grAllThis->Set(1);

  TGraphErrors* grSystThis = dynamic_cast<TGraphErrors*>(grSyst->Clone());
  grSystThis->SetFillColorAlpha(kRed, 0.);
  grSystThis->SetFillStyle(0);
  grSystThis->SetMarkerColor(kRed);
  grSystThis->SetLineColor(kRed);
  grSystThis->Set(1);

  TLine* linePdg = new TLine(pdgValue, 0, pdgValue, dataPoints.size());
  linePdg->SetLineColor(kBlue);
  linePdg->SetLineWidth(3);

  TBox* boxPdg = new TBox(pdgValue - pdgError, 0, pdgValue + pdgError, dataPoints.size());
  boxPdg->SetFillStyle(1000);
  boxPdg->SetFillColorAlpha(kBlue, 0.2);
  boxPdg->SetLineColor(kBlue);

  if(isCombineStatSyst) {
    grAll->GetXaxis()->SetLimits(162, 234);
    grAll->Draw("APE");
    linePdg->Draw("same");
    boxPdg->Draw("same");
    grAll->Draw("PE same");
    grAllThis->Draw("PE same");
  } else {
    grSyst->GetXaxis()->SetLimits(162, 234);
    grSyst->Draw("AP5");
    linePdg->Draw("same");
    boxPdg->Draw("same");
    grSyst->Draw("P5 same");
    grAll->Draw("PE same");
    grSystThis->Draw("P5 same");
    grAllThis->Draw("PE same");
  }

  TLegend* leg = new TLegend(0.405, 0.84, 0.645, 0.97);
  leg->SetTextSize(0.043);
  leg->SetHeader("ALICE Preliminary");
  leg->AddEntry(boxPdg, "PDG average (2025)", "LF");
  leg->Draw("same");

  grAll->GetYaxis()->Set(dataPoints.size(), 0, dataPoints.size());
  grAll->GetYaxis()->SetRangeUser(0, dataPoints.size());
  grSyst->GetYaxis()->Set(dataPoints.size(), 0, dataPoints.size());
  grSyst->GetYaxis()->SetRangeUser(0, dataPoints.size());
  for(int iDp=0, nDps=dataPoints.size(); iDp<nDps; ++iDp) {
    const auto& dp = dataPoints.at(iDp);
    grAll->GetYaxis()->SetBinLabel(grAll->GetYaxis()->FindBin(iDp+0.5), formLabel(dp).c_str());
    grSyst->GetYaxis()->SetBinLabel(grAll->GetYaxis()->FindBin(iDp+0.5), formLabel(dp).c_str());
  }

  grAll->RemovePoint(0);
  grSyst->RemovePoint(0);

  cc->Print("cc.pdf", "pdf");
}
