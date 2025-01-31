void lambdac_minv(const std::string& fileName = "AnalysisResults.merged.issel.root") {
  gROOT->Macro( "style_1Dhisto.cc" );

  TFile* fileIn = TFile::Open(fileName.c_str());

  const int rebinfactor{1};

  TCanvas* cc = new TCanvas("ccEff", "efficiency", 1200, 800);

  TH1F* histoSig = fileIn->Get<TH1F>("hf-task-lc/MC/reconstructed/signal/hMassRecSig");
  TH1F* histoBg = fileIn->Get<TH1F>("hf-task-lc/MC/reconstructed/background/hMassRecBg");
  if(histoSig == nullptr || histoBg == nullptr) {
    throw std::runtime_error("histoSig == nullptr || histoBg == nullptr");
  }
  if(rebinfactor != 1) {
    histoSig->Rebin(rebinfactor);
    histoBg->Rebin(rebinfactor);
  }

  TH1F* histoAll = (TH1F*)histoSig->Clone();
  histoAll->SetName("hAll");
  histoAll->SetTitle("hAll");
  histoAll->GetYaxis()->SetTitle("Entries");
  histoAll->Add(histoBg);
  histoAll->SetLineWidth(2);
  histoAll->SetTitle("");
  histoAll->SetStats(0);
  cc->cd();
  histoAll->Draw();

  cc->Print("ccAll.pdf", "pdf");

  fileIn->Close();
}
