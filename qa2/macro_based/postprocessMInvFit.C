#include "/home/oleksii/alidir/macros_on_git/qa2/exe_based/HelperGeneral.hpp"

using namespace HelperGeneral;

void postprocessMInvFit() {
  const std::string fileName{"2/RawYields_Lc/RawYields_Lc.NPgt0.01.root"};
  const int canvasNumber{4}; // starts from 1
  const bool isResidual{false};

  gROOT->Macro("/home/oleksii/alidir/working/hpApprovals/postprocessMInvFit.style.cc");
  const double titleSize{0.05};
  const double legendSize{0.04};

  TFile* fileIn = OpenFileWithNullptrCheck(fileName.c_str());

  std::string canvasName, histoName, mainFuncName;
  if(!isResidual) {
    canvasName = "canvasMass0";
    histoName = "data_c";
    mainFuncName = "Tot_c";
  } else {
    canvasName = "canvasResiduals0";
    histoName = "resid_data_c_Bkg_c";
    mainFuncName = "mSgnPdf_Norm[mass]";
  }

  TCanvas* cc = GetObjectWithNullptrCheck<TCanvas>(fileIn, canvasName);

  TPad* pad = (TPad*)cc->GetPad(canvasNumber);
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
  histo->GetXaxis()->SetTitle("#it{M}(pK#pi) (GeV/#it{c}^{2})");
  histo->GetXaxis()->SetTitleSize(titleSize);
  histo->GetXaxis()->SetLabelSize(titleSize);
  histo->GetXaxis()->SetTitleOffset(1.1);
  histo->GetYaxis()->SetTitle(("Counts per " + to_string_with_precision(binWidth, 0) + " MeV/#kern[-0.0]{#it{c}}^{2}").c_str());
  histo->GetYaxis()->SetTitleSize(titleSize);
  histo->GetYaxis()->SetLabelSize(titleSize);
  histo->GetYaxis()->SetTitleOffset(1.1);
  histo->GetYaxis()->SetNdivisions(206);
  histo->GetYaxis()->SetMaxDigits(3);
  if(!isResidual) {
    histo->SetMinimum(0.);
    histo->SetMaximum(1.12 * histo->GetYaxis()->GetXmax()); // FIXME ad. hoc. to avoid overlap with left text
  }
  histo->SetMarkerSize(1.4);

  RooCurve* fitBkg = !isResidual ? dynamic_cast<RooCurve*>(list1->FindObject("Bkg_c")) : nullptr;
  if(fitBkg == nullptr && !isResidual) throw std::runtime_error("fitBkg == nullptr");

  RooCurve* fitTotal = dynamic_cast<RooCurve*>(list1->FindObject(mainFuncName.c_str()));
  if(fitTotal == nullptr) throw std::runtime_error("fitTotal == nullptr");

  RooCurve* fitSignal = !isResidual ? dynamic_cast<RooCurve*>(list1->FindObject("mSgnPdf_Norm[mass]")) : nullptr;
  if(fitSignal == nullptr && !isResidual) throw std::runtime_error("fitSignal == nullptr");
  if(!isResidual) {
    fitSignal->SetFillColorAlpha(kAzure+9, 0.4);
    fitSignal->SetLineStyle(0);
    fitSignal->SetLineWidth(0);
    fitSignal->SetLineColor(0);
  }

  TH1* hSignalYield = GetObjectWithNullptrCheck<TH1>(fileIn, "hRawYieldsSignal");
  const double signal = hSignalYield->GetBinContent(canvasNumber);
  const double signalError = hSignalYield->GetBinError(canvasNumber);
  const double decayTimeLo = hSignalYield->GetBinLowEdge(canvasNumber);
  const double decayTimeUp = hSignalYield->GetBinLowEdge(canvasNumber + 1);

  TPaveText* leftText = new TPaveText(0.15, 0.64, 0.25, 0.92, "brNDC");
  leftText->SetTextAlign(12);
  leftText->SetFillColor(0);
  leftText->SetTextSize(legendSize);
  leftText->SetTextFont(42);
  leftText->AddText("ALICE Preliminary");
  leftText->AddText("pp,#kern[0.3]{#sqrt{#it{s}}} = 13.6 TeV");
  leftText->AddText("#it{L}_{int} = 4.9 pb^{-1}");
  leftText->AddText("#Lambda_{c}^{+}#kern[0.5]{#rightarrow} pK^{-}#pi^{+} + c.c.");
  leftText->AddText((to_string_with_precision(decayTimeLo, 1) + " <#kern[1.0]{#it{t}} < " + to_string_with_precision(decayTimeUp, 1) + " ps").c_str());
  leftText->AddText("3 <#kern[0.3]{#it{p}_{T}} < 20 GeV/#it{c}");

  TLegend* legend = new TLegend(0.15, 0.28, 0.35, 0.45);
  legend->SetBorderSize(0);
  legend->SetTextSize(legendSize);
  legend->SetTextFont(42);
  if(!isResidual) {
    legend->AddEntry(histo, "Data", "PE");
    legend->AddEntry(fitTotal, "Total fit function", "L");
    legend->AddEntry(fitBkg, "Combinatorial bkg", "L");
    legend->AddEntry(fitSignal, "Signal", "F");
  } else {
    legend->AddEntry(histo, "residuals", "PE");
    legend->AddEntry(fitTotal, "signal", "L");
  }

  TCanvas* c1 = new TCanvas("c1", "", 1000, 1000);
  c1->SetTicks(1, 1);
  histo->Draw("AP");
  if(!isResidual) {
    fitBkg->Draw("same");
    fitSignal->Draw("F same");
  }
  fitTotal->Draw("same");
  leftText->Draw("same");
  legend->Draw("same");

  c1->Print((canvasName + ".eps").c_str(), "eps");
}
