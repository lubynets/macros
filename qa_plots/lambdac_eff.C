void lambdac_eff(const std::string& fileName = "AnalysisResults.merged.nosel.root") {
  gROOT->Macro( "style_1Dhisto.cc" );

  TFile* fileIn = TFile::Open(fileName.c_str());

  struct SignalSpecies{
    std::string dirname_;
    std::string histonameend_;
    Color_t color_;
  };

  std::vector<SignalSpecies> signal_species {
    {"signal",    ""         , kBlue    },
    {"prompt",    "Prompt"   , kRed     },
    {"nonprompt", "NonPrompt", kGreen+2 }
  };

  TFile* fileOut = TFile::Open("fileOut.root", "recreate");

  const int rebinfactor{1};

  TCanvas* cc = new TCanvas("ccEff", "efficiency", 1200, 800);
  TLegend* leg = new TLegend(0.2, 0.7, 0.4, 0.82);
  leg->SetBorderSize(0);

  bool isFirstHisto{true};
  for(auto& ss : signal_species) {
    TH1F* histoRec = fileIn->Get<TH1F>(("hf-task-lc/MC/reconstructed/" + ss.dirname_ + "/hPtRecSig" + ss.histonameend_).c_str());
    TH1F* histoGen = fileIn->Get<TH1F>(("hf-task-lc/MC/generated/" + ss.dirname_ + "/hPtGen" + ss.histonameend_).c_str());
    if(histoRec == nullptr || histoGen == nullptr) {
      throw std::runtime_error("histoRec == nullptr || histoGen == nullptr");
    }
    if(rebinfactor != 1) {
      histoRec->Rebin(rebinfactor);
      histoGen->Rebin(rebinfactor);
    }
    TH1F* histoEff = (TH1F*)histoRec->Clone();
    histoEff->SetName(("hEffSig" + ss.histonameend_).c_str());
    histoEff->SetTitle(ss.dirname_.c_str());
    histoEff->GetYaxis()->SetTitle("#varepsilon, %");
    histoEff->Sumw2();
    histoEff->Divide(histoGen);
    histoEff->Scale(100);
    histoEff->GetXaxis()->SetRangeUser(0, 15);
    histoEff->GetYaxis()->SetRangeUser(0, 20);
    histoEff->Write();
    histoEff->SetLineColor(ss.color_);
    histoEff->SetLineWidth(2);
    histoEff->SetTitle("");
    histoEff->SetStats(0);
    if(isFirstHisto) histoEff->Draw("");
    else             histoEff->Draw("same");
    leg->AddEntry(histoEff, ss.dirname_.c_str(), "L");
    isFirstHisto = false;
  }

  leg->Draw();
  cc->Write();
  cc->Print("ccEff.pdf");
  fileOut->Close();
  fileIn->Close();
}
