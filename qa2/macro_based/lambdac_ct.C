void lambdac_ct(const std::string& fileName = "AnalysisResults.merged.issel.root") {
  gROOT->Macro( "style_1Dhisto.cc" );

  TFile* fileIn = TFile::Open(fileName.c_str());

  struct SignalSpecies{
    std::string dirname_;
    std::string histonameend_;
    Color_t color_;
  };

  std::vector<SignalSpecies> signal_species {
//     {"signal"    , "Sig"         , kBlue    },
    {"prompt"    , "SigPrompt"   , kRed     },
    {"nonprompt" , "SigNonPrompt", kGreen+2 },
    {"background", "Bg"          , kBlack   }
  };

  TFile* fileOut = TFile::Open("fileOut.root", "recreate");

  const int rebinfactor{1};

  TCanvas* cc = new TCanvas("ccCt", "ct", 1200, 800);
  cc->SetLogy();
  TLegend* leg = new TLegend(0.6, 0.7, 0.8, 0.82);
  leg->SetBorderSize(0);

  bool isFirstHisto{true};
  double integral{0};
  double scalefactor{1};
  double ymax{0};
  for(auto& ss : signal_species) {
    TH1F* histoCt = fileIn->Get<TH1F>(("hf-task-lc/MC/reconstructed/" + ss.dirname_ + "/hCtRec" + ss.histonameend_).c_str());
    if(histoCt == nullptr) {
      throw std::runtime_error("histoCt == nullptr");
    }
    if(rebinfactor != 1) {
      histoCt->Rebin(rebinfactor);
    }
    if(isFirstHisto) integral = histoCt->GetEntries();
    else             scalefactor = integral / histoCt->GetEntries();
    histoCt->SetName(("hCt" + ss.histonameend_).c_str());
    histoCt->SetTitle(ss.dirname_.c_str());
    histoCt->GetYaxis()->SetTitle("a. u.");
    histoCt->Sumw2();
    if(!isFirstHisto) histoCt->Scale(scalefactor);
    histoCt->Write();
    histoCt->SetLineColor(ss.color_);
    histoCt->SetLineWidth(2);
    histoCt->SetTitle("");
    histoCt->SetStats(0);
    if(isFirstHisto) histoCt->Draw("");
    else             histoCt->Draw("same");
    leg->AddEntry(histoCt, ss.dirname_.c_str(), "L");
    isFirstHisto = false;
  }

  leg->Draw();
  cc->Write();
  cc->Print("ccCt.pdf");
  fileOut->Close();
  fileIn->Close();
}
