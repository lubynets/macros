#include "/home/oleksii/alidir/macros_on_git/qa2/exe_based/HelperGeneral.hpp"

using namespace HelperGeneral;

void postprocessFit() {
  const std::string fileName{"2/RawYields_Lc/RawYields_Lc.NPgt0.01.root"};
  const int canvasNumber{4}; // starts from 1
  const bool isResidual{false};


  gROOT->Macro("/home/oleksii/alidir/macros_on_git/qa2/macro_based/styles/postprocessFit.style.cc");
  const double titleSize{0.05};
  const double legendSize{0.04};

  TFile* fileIn = OpenFileWithNullptrCheck(fileName.c_str());

  std::string canvasName, histoName, sgnFuncName;
  if(!isResidual) {
    canvasName = "canvasMass0";
    histoName = "data_c";
    sgnFuncName = "Tot_c";
  } else {
    canvasName = "canvasResiduals0";
    histoName = "resid_data_c_Bkg_c";
    sgnFuncName = "mSgnPdf_Norm[mass]";
  }

  TCanvas* cc = GetObjectWithNullptrCheck<TCanvas>(fileIn, canvasName);

  TPad *pad = (TPad*)cc->GetPad(canvasNumber);
  if(pad == nullptr) throw std::runtime_error("pad == nullptr");
  pad->cd();
  pad->Modified();
  pad->Update();

  TList* list1 = pad->GetListOfPrimitives();
  if(list1 == nullptr) throw std::runtime_error("list1 == nullptr");

  RooHist* histo = dynamic_cast<RooHist*>(list1->FindObject(histoName.c_str()));
  if(histo == nullptr) throw std::runtime_error("histo == nullptr");
  const double binWidth = (histo->GetX()[2] - histo->GetX()[1]) * 1000;
  histo->SetTitle("");
  histo->GetXaxis()->SetRangeUser(2.12, 2.42);
  for (int i = 0; i < histo->GetN(); i++) {
    histo->SetPointEXlow(i, 0);
    histo->SetPointEXhigh(i, 0);
  }
  histo->GetXaxis()->SetTitle("M_{pK#pi} (GeV/#it{c}^{2})");
  histo->GetXaxis()->SetTitleSize(titleSize);
  histo->GetXaxis()->SetLabelSize(titleSize);
  histo->GetXaxis()->SetTitleOffset(1.1);
  histo->GetYaxis()->SetTitle(("Counts per " + to_string_with_precision(binWidth, 0) + " MeV/#it{c}^{2}").c_str());
  histo->GetYaxis()->SetTitleSize(titleSize);
  histo->GetYaxis()->SetLabelSize(titleSize);
  histo->GetYaxis()->SetTitleOffset(2.1);
  histo->GetYaxis()->SetNdivisions(206);
  histo->SetMarkerSize(1.4);

  RooCurve* fitBkg = !isResidual ? dynamic_cast<RooCurve*>(list1->FindObject("Bkg_c")) : nullptr;
  if(fitBkg == nullptr && !isResidual) throw std::runtime_error("fitBkg == nullptr");

  RooCurve* fitTotal = dynamic_cast<RooCurve*>(list1->FindObject(sgnFuncName.c_str()));
  if(fitTotal == nullptr) throw std::runtime_error("fitTotal == nullptr");

  TH1* hSignalYield = GetObjectWithNullptrCheck<TH1>(fileIn, "hRawYieldsSignal");
  const double signal = hSignalYield->GetBinContent(canvasNumber);
  const double signalError = hSignalYield->GetBinError(canvasNumber);
  const double decayTimeLo = hSignalYield->GetBinLowEdge(canvasNumber);
  const double decayTimeUp = hSignalYield->GetBinLowEdge(canvasNumber + 1);

  TH1* hChi2Ndf = GetObjectWithNullptrCheck<TH1>(fileIn, "hRawYieldsChiSquareTotal");
  const double chi2ndf = hChi2Ndf->GetBinContent(canvasNumber);

  TPaveText* leftText = new TPaveText(0.35, 0.82, 0.45, 0.89, "brNDC");
  leftText->SetFillColor(0);
  leftText->SetTextSize(legendSize);
  leftText->SetTextFont(62);
  leftText->AddText("#Lambda_{c}^{+}#rightarrow pK^{-}#pi^{+} + c.c.");
  leftText->AddText("pp #sqrt{s} = 13.6 TeV");

  TPaveText* rightText = new TPaveText(0.70, 0.82, 0.90, 0.89, "brNDC");
  rightText->SetFillColor(0);
  rightText->SetTextSize(legendSize);
  rightText->SetTextFont(62);
  rightText->AddText(("S = " + to_string_with_precision(signal, 0) + "#pm " + to_string_with_precision(signalError, 0)).c_str());
  rightText->AddText(("#chi^{2}/ndf = " + to_string_with_precision(chi2ndf, 2)).c_str());

  if(!isResidual) histo->SetTitle(("t#in (" + to_string_with_precision(decayTimeLo, 1) + "#; " + to_string_with_precision(decayTimeUp, 1) + ") (ps)").c_str());

  TLegend* legend = new TLegend(0.25, 0.62, 0.45, 0.75);
  legend->SetBorderSize(0);
  legend->SetTextSize(legendSize);
  legend->SetTextFont(62);
  if(!isResidual) {
    legend->AddEntry(histo, "Data", "PE");
    legend->AddEntry(fitTotal, "Fit SIG + BG", "L");
    legend->AddEntry(fitBkg, "Fit BG", "L");
  } else {
    legend->AddEntry(histo, "Data #minus Fit BG", "PE");
    legend->AddEntry(fitTotal, "Fit SIG", "L");
  }

  TCanvas* c1 = new TCanvas("c1", "", 1000, 1000);
  c1->SetTicks(1, 1);
  histo->Draw("AP");
  if(!isResidual) fitBkg->Draw("same");
  fitTotal->Draw("same");
  if(!isResidual) {
    leftText->Draw("same");
    rightText->Draw("same");
  }
  legend->Draw("same");

  c1->Print((canvasName + ".pdf").c_str(), "pdf");
}
