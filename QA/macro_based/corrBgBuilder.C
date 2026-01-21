template<typename T>
inline std::string to_string_with_precision(const T a_value, const int n=2) {
  std::ostringstream out;
  out.precision(n);
  out << std::fixed << a_value;
  return out.str();
}

void corrBgBuilder() {
  gStyle->SetOptStat(0);
  const std::string fileName = "/home/oleksii/alidir/working/correlBG/LHC25e4/corrBg_qa.HL.mc.HF_LHC25e4_All.579594.ctbin3.root";
  TFile* fileIn = TFile::Open(fileName.c_str(), "read");

//   const std::vector<float> ptRanges = {3,/* 2, 3, 4, 5, 8, 12, */20};
  const std::vector<float> ptRanges = {3, 4};
  const std::vector<float> ctRanges = {0.0, 20.0}; const int ctPrecision = 1;
  const int rebinFactor{4};

  for(int iPt=0, nPts=ptRanges.size()-1; iPt<nPts; ++iPt) {
    for(int iCt=0, nCts=ctRanges.size()-1; iCt<nCts; ++iCt) {
      const std::string hPath = "pT_" + to_string_with_precision(ptRanges.at(iPt), 0) + "_" + to_string_with_precision(ptRanges.at(iPt+1), 0) + "/T_" + to_string_with_precision(ctRanges.at(iCt), ctPrecision) + "_" + to_string_with_precision(ctRanges.at(iCt+1), ctPrecision) + "/";
      TH1* hSig = fileIn->Get<TH1>((hPath + "hMass_LcToPKPiSig").c_str());
      TH1* hCorrBg = fileIn->Get<TH1>((hPath + "hMass_bkgSum").c_str());
      if(hSig == nullptr || hCorrBg == nullptr) throw std::runtime_error("hSig == nullptr || hCorrBg == nullptr for " + hPath);

      if(rebinFactor != 1) {
        hSig->Rebin(rebinFactor);
        hCorrBg->Rebin(rebinFactor);
      }
      hSig->SetLineColor(kOrange+2);
      hCorrBg->SetLineColor(kGray+2);
      hSig->SetLineWidth(3);
      hCorrBg->SetLineWidth(3);

      TCanvas ccSignal("ccSignal", "");
      TCanvas ccCorrBg("ccCorrBg", "");
      TCanvas ccBoth("ccBoth", "");
      for(const auto& cc : {&ccSignal, &ccCorrBg, &ccBoth}) {
        cc->SetCanvasSize(1000, 1000);
//         cc->SetLogy();
      }

      TPaveText legText(0.18, 0.78, 0.35, 0.88, "brNDC");
      legText.AddText(("#it{p}_{T} #in (" + to_string_with_precision(ptRanges.at(iPt), 0) + "; " + to_string_with_precision(ptRanges.at(iPt+1), 0) + ") GeV/#it{c}").c_str());
      legText.AddText(("t #in (" + to_string_with_precision(ctRanges.at(iCt)) + "; " + to_string_with_precision(ctRanges.at(iCt+1)) + ") ps").c_str());
      legText.SetFillColor(0);
      legText.SetTextSize(0.03);
      legText.SetTextFont(62);

      ccSignal.cd();
      hSig->Draw();
      legText.Draw("same");

      ccCorrBg.cd();
      hCorrBg->Draw();
      legText.Draw("same");

      ccBoth.cd();
      hSig->Draw();
      hCorrBg->Draw("same");
      legText.Draw("same");

      for(const auto& cc : {&ccSignal, &ccCorrBg, &ccBoth}) {
        cc->Print((static_cast<std::string>(cc->GetName()) + "_pt" + std::to_string(iPt) + "_ct" + std::to_string(iCt) + ".pdf").c_str(), "pdf");
      }
    }
  }
}
