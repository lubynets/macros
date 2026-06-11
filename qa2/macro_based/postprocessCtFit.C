#include "/home/oleksii/alidir/macros_on_git/qa2/exe_based/HelperGeneral.hpp"

using namespace HelperGeneral;

void drawFrame(TH1* histoIn, bool isLogy=false, double filledAreaX = 0.8, double filledAreaY = 0.7);

void postprocessCtFit() {
  const double legendSize{0.045};
  gROOT->Macro("/home/oleksii/alidir/working/hpApprovals/postprocessCtFit.style.cc");

  const std::string fileName{"/home/oleksii/alidir/working/minvfitDebug/1stpol3/median/NP/W/ctfit.root"};

  TFile* fileIn = OpenFileWithNullptrCheck(fileName.c_str());

  TCanvas* ccFit = GetObjectWithNullptrCheck<TCanvas>(fileIn, "ccFit");
  TList* list1 = ccFit->GetListOfPrimitives();
  if(list1 == nullptr) throw std::runtime_error("list1 == nullptr");

  TH1* histoYield = dynamic_cast<TH1*>(list1->FindObject("hCorrYieldsPrompt"));
  if(histoYield == nullptr) throw std::runtime_error("histoYield == nullptr");
  histoYield->GetYaxis()->SetTitleOffset(1.0);
  histoYield->GetXaxis()->SetDecimals(true);

  TF1* fitFunc = dynamic_cast<TF1*>(list1->FindObject("fitFunc"));
  if(fitFunc == nullptr) throw std::runtime_error("fitFunc == nullptr");
  fitFunc->SetLineColor(kBlue);
  fitFunc->SetLineStyle(7);

  TCanvas* ccFitPP = new TCanvas("ccFitPP", "");
  ccFitPP->SetCanvasSize(ccFit->GetWw(), ccFit->GetWh());
  ccFitPP->SetLogy();
  ccFitPP->SetTicks(1, 1);

  drawFrame(histoYield, true);
  histoYield->Draw("same");
  fitFunc->Draw("same");

  TPaveText* rightText = new TPaveText(0.60, 0.64, 0.75, 0.92, "brNDC");
  rightText->SetTextAlign(12);
  rightText->SetFillColor(0);
  rightText->SetTextSize(legendSize);
  rightText->SetTextFont(42);
  rightText->AddText("ALICE Preliminary");
  rightText->AddText("pp,#kern[0.3]{#sqrt{#it{s}}} = 13.6 TeV");
  rightText->AddText("#it{L}_{int} = 4.9 pb^{-1}");
  rightText->AddText("#Lambda_{c}^{+}#kern[0.5]{#rightarrow} pK^{-}#pi^{+} + c.c.");
  rightText->AddText("3 <#kern[0.3]{#it{p}_{T}} < 20 GeV/#it{c}");
  rightText->Draw("same");

  TPaveText* tauText = new TPaveText(0.60, 0.55, 0.75, 0.60, "brNDC");
  tauText->SetTextAlign(12);
  tauText->SetFillColor(0);
  tauText->SetTextSize(legendSize);
  tauText->SetTextFont(42);
  tauText->AddText(("#tau_{#Lambda_{c}} = " +
                                            to_string_with_precision(fitFunc->GetParameter(1)*1000, 1) +
                                            "#kern[0.7]{#pm} " +
                                            to_string_with_precision(fitFunc->GetParError(1)*1000, 1) +
                                            "(stat.) fs").c_str());
  tauText->Draw("same");

  TLegend* leg = new TLegend(0.2, 0.3, 0.4, 0.5);
  leg->SetBorderSize(0);
  leg->SetTextSize(legendSize);
  leg->SetTextFont(42);
  leg->AddEntry(histoYield, "Data", "PE");
  leg->AddEntry(fitFunc, "Fit", "L");
  leg->Draw("same");

  ccFitPP->Print("ccFit.pdf(", "pdf");

  // ========================================================================================================

  TCanvas* ccRatio = GetObjectWithNullptrCheck<TCanvas>(fileIn, "ccRatio");
  list1 = ccRatio->GetListOfPrimitives();
  if(list1 == nullptr) throw std::runtime_error("list1 == nullptr");

  TH1* histoRatio = dynamic_cast<TH1*>(list1->FindObject("hCorrYieldsPrompt"));
  if(histoRatio == nullptr) throw std::runtime_error("histoRatio == nullptr");
  histoRatio->GetYaxis()->SetTitleOffset(0.62);
  histoRatio->GetYaxis()->SetRangeUser(0.9401, 1.045);
  histoRatio->GetXaxis()->SetTickLength(2*histoYield->GetXaxis()->GetTickLength());

  TF1* oneLine = dynamic_cast<TF1*>(list1->FindObject("oneline"));
  if(oneLine == nullptr) throw std::runtime_error("oneLine == nullptr");

  TCanvas* ccRatioPP = new TCanvas("ccRatioPP", "");
  ccRatioPP->SetCanvasSize(ccRatio->GetWw(), ccRatio->GetWh());
  ccRatioPP->SetBottomMargin(2*ccFitPP->GetBottomMargin());
  ccRatioPP->SetTopMargin(2*ccFitPP->GetTopMargin());
  ccRatioPP->SetTicks(1, 1);

  drawFrame(histoRatio);
  histoRatio->Draw("same");
  oneLine->Draw("same");

  ccRatioPP->Print("ccFit.pdf", "pdf");

  // ========================================================================================================

  TCanvas* ccResidual = GetObjectWithNullptrCheck<TCanvas>(fileIn, "ccResidual");
  list1 = ccResidual->GetListOfPrimitives();
  if(list1 == nullptr) throw std::runtime_error("list1 == nullptr");

  TH1* histoResidual = dynamic_cast<TH1*>(list1->FindObject("hCorrYieldsPrompt"));
  if(histoResidual == nullptr) throw std::runtime_error("histoResidual == nullptr");
  histoResidual->GetYaxis()->SetTitleOffset(0.62);
  histoResidual->GetXaxis()->SetTickLength(2*histoYield->GetXaxis()->GetTickLength());

  TF1* zeroLine = dynamic_cast<TF1*>(list1->FindObject("oneline"));
  if(zeroLine == nullptr) throw std::runtime_error("zeroLine == nullptr");

  TCanvas* ccResidualPP = new TCanvas("ccResidualPP", "");
  ccResidualPP->SetCanvasSize(ccResidual->GetWw(), ccResidual->GetWh());
  ccResidualPP->SetBottomMargin(2*ccFitPP->GetBottomMargin());
  ccResidualPP->SetTopMargin(2*ccFitPP->GetTopMargin());
  ccResidualPP->SetTicks(1, 1);

  drawFrame(histoResidual);
  histoResidual->Draw("same");
  zeroLine->Draw("same");

  ccResidualPP->Print("ccFit.pdf)", "pdf");

}

void drawFrame(TH1* histoIn, bool isLogy, const double filledAreaX, const double filledAreaY) {

  const double xmin = histoIn->GetXaxis()->GetXmin();
  const double xmax = histoIn->GetXaxis()->GetXmax();
  double ymin = histoIn->GetMinimum();
  double ymax = histoIn->GetMaximum();
  if(isLogy) {
    ymin = std::log10(ymin);
    ymax = std::log10(ymax);
  }

  const double xdiff = xmax - xmin;
  const double ydiff = ymax - ymin;

  const double xlo = xmin - (1-filledAreaX) * xdiff / 2;
  const double xhi = xmax + (1-filledAreaX) * xdiff / 2;
  double ylo = ymin - (1-filledAreaY) * ydiff / 2;
  double yhi = ymax + (1-filledAreaY) * ydiff / 2;
  if(isLogy) {
    ylo = std::pow(10, ylo);
    yhi = std::pow(10, yhi);
  }

  TH1* hFrame = gPad->DrawFrame(xlo, ylo, xhi, yhi);
  hFrame->GetXaxis()->SetTitle(histoIn->GetXaxis()->GetTitle());
  hFrame->GetYaxis()->SetTitle(histoIn->GetYaxis()->GetTitle());
  hFrame->GetXaxis()->SetTitleSize(histoIn->GetXaxis()->GetTitleSize());
  hFrame->GetYaxis()->SetTitleSize(histoIn->GetYaxis()->GetTitleSize());
  hFrame->GetXaxis()->SetLabelSize(histoIn->GetXaxis()->GetLabelSize());
  hFrame->GetYaxis()->SetLabelSize(histoIn->GetYaxis()->GetLabelSize());
  hFrame->GetXaxis()->SetTitleOffset(histoIn->GetXaxis()->GetTitleOffset());
  hFrame->GetYaxis()->SetTitleOffset(histoIn->GetYaxis()->GetTitleOffset());
  hFrame->GetXaxis()->SetLabelOffset(histoIn->GetXaxis()->GetLabelOffset());
  hFrame->GetYaxis()->SetLabelOffset(histoIn->GetYaxis()->GetLabelOffset());
  hFrame->GetXaxis()->SetTickLength(histoIn->GetXaxis()->GetTickLength());
  hFrame->GetXaxis()->SetDecimals();
  hFrame->GetYaxis()->SetDecimals(!isLogy);
}
