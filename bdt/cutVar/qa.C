TH1* evalTH1DiffNorm(const TH1* h1, const TH1* h2); // bin-by-bin difference of mean values normalized by Barlow-calculated sigma
TH1* evalTH1Sigma(const TH1* h);
TH1* evalTH1SigmaDiff(const TH1* h1, const TH1* h2); // difference of errors in quadrature

void qa() {
  TFile* fileIs = TFile::Open("/home/oleksii/alidir/working/cutVar/testZeroEff/isStatErr/CutVarLc.merged.root", "read");
  TFile* fileNo = TFile::Open("/home/oleksii/alidir/working/cutVar/testZeroEff/noStatErr/CutVarLc.merged.root", "read");
  TFile* fileSys = TFile::Open("/home/oleksii/alidir/working/cutVar/testZeroEff/hSyst.root", "read");

  TH1* histoChi2Is = fileIs->Get<TH1>("hChi2OverNdf");
  TH1* histoChi2No = fileNo->Get<TH1>("hChi2OverNdf");
  TH1* histoSyst = fileSys->Get<TH1>("hCorrYieldsPrompt");
  histoSyst->SetLineColor(kBlue);

  histoChi2Is->SetLineColor(kRed);
  histoChi2No->SetLineColor(kBlue);
  histoChi2No->GetYaxis()->SetRangeUser(0., histoChi2No->GetMaximum()*1.2);

  TLegend legChi2(0.6, 0.7, 0.85, 0.85);
  legChi2.SetBorderSize(0);
  legChi2.AddEntry(histoChi2Is, "with eff. stat. err.", "L");
  legChi2.AddEntry(histoChi2No, "w/o eff. stat. err.", "L");

  TCanvas cChi2("cChi2", "", 1200, 800);
  histoChi2No->Draw("");
  histoChi2Is->Draw("same");
  legChi2.Draw("same");
  cChi2.Print("chi2.pdf", "pdf");

  TH1* histoYieldPromptIs = fileIs->Get<TH1>("hCorrYieldsPrompt");
  TH1* histoYieldPromptNo = fileNo->Get<TH1>("hCorrYieldsPrompt");

  TH1* hDiffNorm = evalTH1DiffNorm(histoYieldPromptIs, histoYieldPromptNo);
  hDiffNorm->GetYaxis()->SetTitle("(Y_{1}-Y_{2}) / #sqrt{#sigma_{Y1}^{2} - #sigma_{Y2}^{2}}");
  TCanvas cYieldDiffNorm("cYieldDiffNorm", "", 1200, 800);
  hDiffNorm->Draw("");
  cYieldDiffNorm.Print("yieldDiffNorm.pdf", "pdf");

  TH1* hErrIs = evalTH1Sigma(histoYieldPromptIs);
  TH1* hErrNo = evalTH1Sigma(histoYieldPromptNo);
  hErrIs->SetLineColor(kRed);
  hErrNo->SetLineColor(kBlue);
  hErrIs->GetYaxis()->SetTitle("#sigma_{Y}");

  TCanvas cSigma("cSigma", "", 1200, 800);
  hErrIs->Draw("");
  hErrNo->Draw("same");
  cSigma.Print("sigma.pdf", "pdf");

  TH1* hSigmaDiff = evalTH1SigmaDiff(histoYieldPromptIs, histoYieldPromptNo);
  hSigmaDiff->GetYaxis()->SetTitle("");
  hSigmaDiff->SetLineColor(kRed);
  hSigmaDiff->GetYaxis()->SetRangeUser(0., 1.2*std::max(hSigmaDiff->GetMaximum(), histoSyst->GetMaximum()));

  TLegend legSigmaDiff(0.6, 0.7, 0.85, 0.85);
  legSigmaDiff.SetBorderSize(0);
  legSigmaDiff.AddEntry(hSigmaDiff, "#sqrt{#sigma_{Y}^{2}(with stat) - #sigma_{Y}^{2}(w/o stat)}", "L");
  legSigmaDiff.AddEntry(histoSyst, "#sigma_{Y} syst.", "L");

  TCanvas cSigmaDiff("cSigmaDiff", "", 1200, 800);
  hSigmaDiff->Draw("");
  histoSyst->Draw("same");
  legSigmaDiff.Draw("same");
  cSigmaDiff.Print("sigmaDiff.pdf", "pdf");

  fileSys->Close();
  fileNo->Close();
  fileIs->Close();
}

TH1* evalTH1DiffNorm(const TH1* h1, const TH1* h2) {
  TH1* hResult = dynamic_cast<TH1*>(h1->Clone());
  hResult->Reset();

  for(int iBin=1, nBins=h1->GetNbinsX(); iBin<=nBins; ++iBin) {
    const double v1 = h1->GetBinContent(iBin);
    const double v2 = h2->GetBinContent(iBin);
    const double e1 = h1->GetBinError(iBin);
    const double e2 = h2->GetBinError(iBin);
    double diff = v1 - v2;
    const double err = std::sqrt(std::fabs(e1*e1 - e2*e2));
    diff /= err;
    hResult->SetBinContent(iBin, diff);
  }

  return hResult;
}

TH1* evalTH1Sigma(const TH1* h) {
  TH1* hResult = dynamic_cast<TH1*>(h->Clone());
  hResult->Reset();

  for(int iBin=1, nBins=h->GetNbinsX(); iBin<=nBins; ++iBin) {
    hResult->SetBinContent(iBin, h->GetBinError(iBin));
  }

  return hResult;
}

TH1* evalTH1SigmaDiff(const TH1* h1, const TH1* h2) {
  TH1* hResult = dynamic_cast<TH1*>(h1->Clone());
  hResult->Reset();

  for(int iBin=1, nBins=h1->GetNbinsX(); iBin<=nBins; ++iBin) {
    const double e1 = h1->GetBinError(iBin);
    const double e2 = h2->GetBinError(iBin);
    const double diff = std::sqrt(std::fabs(e1*e1 - e2*e2));
    hResult->SetBinContent(iBin, diff);
  }

  return hResult;
}
