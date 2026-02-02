template<typename T>
inline std::string to_string_with_precision(const T a_value, const int n=2) {
  std::ostringstream out;
  out.precision(n);
  out << std::fixed << a_value;
  return out.str();
}

const std::vector<Color_t> palette{kBlue, kRed, kGreen+2, kMagenta, kBlack, kOrange+2, kCyan+2, kViolet-3, kOrange+10};

void corrBgBuilder() {
  gStyle->SetOptStat(0);
  const std::string fileName = "/home/oleksii/alidir/working/correlBG/LHC25e4/corrBg_qa.HL.mc.HF_LHC25e4_All.595825.ctbin3.root";
  TFile* fileIn = TFile::Open(fileName.c_str(), "read");

  const std::vector<float> ptRanges = {1, 2, 3, 4, 5, 8, 12, 20};
//   const std::vector<float> ptRanges = {4, 20};
//   const std::vector<float> ctRanges = {0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.4, 3.6, 5.0}; const int ctPrecision = 1;
  const std::vector<float> ctRanges = {0, 20}; const int ctPrecision = 1;
  const int rebinFactor{4};
  const int smoothFactor{0};
  const std::vector<float> npLowerCutsForCorrBg{/*0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7*/0.01};
  const float npLowerCutsForBoth{0};

  TFile* fileOut = TFile::Open("ccBoth.root", "recreate");

  for(int iPt=0, nPts=ptRanges.size()-1; iPt<nPts; ++iPt) {
    TLegend leg(0.12, 0.55, 0.3, 0.75);
    leg.SetBorderSize(0);

    for(int iCt=0, nCts=ctRanges.size()-1; iCt<nCts; ++iCt) {
      TCanvas ccCorrBg("ccCorrBg", "");
      TCanvas ccBoth("ccBoth", "");
      for(const auto& cc : {&ccCorrBg, &ccBoth}) {
        cc->SetCanvasSize(1000, 1000);
//         cc->SetLogy();
      }

      TPaveText legTextCommon(0.18, 0.74, 0.35, 0.88, "brNDC");
      legTextCommon.AddText(("#it{p}_{T} #in (" + to_string_with_precision(ptRanges.at(iPt), 0) + "; " + to_string_with_precision(ptRanges.at(iPt+1), 0) + ") GeV/#it{c}").c_str());
      legTextCommon.AddText(("t #in (" + to_string_with_precision(ctRanges.at(iCt)) + "; " + to_string_with_precision(ctRanges.at(iCt+1)) + ") ps").c_str());
      legTextCommon.AddText(("smooth = " + std::to_string(smoothFactor)).c_str());
      legTextCommon.SetFillColor(0);
      legTextCommon.SetTextSize(0.03);
      legTextCommon.SetTextFont(62);

      const std::string hPath = "pT_" + to_string_with_precision(ptRanges.at(iPt), 0) + "_" + to_string_with_precision(ptRanges.at(iPt+1), 0) + "/T_" + to_string_with_precision(ctRanges.at(iCt), ctPrecision) + "_" + to_string_with_precision(ctRanges.at(iCt+1), ctPrecision) + "/";

      TH1* hSig = fileIn->Get<TH1>((hPath + "hMass_LcToPKPiSig_NPgt" + to_string_with_precision(npLowerCutsForBoth, 2)).c_str());
      TH1* hCorrBg = fileIn->Get<TH1>((hPath + "hMass_bkgSum_NPgt" + to_string_with_precision(npLowerCutsForBoth, 2)).c_str());
      if(hSig == nullptr || hCorrBg == nullptr) throw std::runtime_error("hSig == nullptr || hCorrBg == nullptr for " + hPath);

      if(rebinFactor != 1) {
        hSig->Rebin(rebinFactor);
        hCorrBg->Rebin(rebinFactor);
      }
      if(smoothFactor != 0) hCorrBg->Smooth(smoothFactor);

      hSig->SetLineColor(kOrange+2);
      hCorrBg->SetLineColor(kGray+2);
      hSig->SetLineWidth(3);
      hCorrBg->SetLineWidth(3);

      ccBoth.cd();
      hSig->Draw();
      hCorrBg->Draw("same");
      TPaveText legTextBoth = legTextCommon;
      legTextBoth.AddText(("bdt(NP) > " + to_string_with_precision(npLowerCutsForBoth, 2)).c_str());
      legTextBoth.Draw("same");

      const std::string ccNameSuffix = "_pt" + std::to_string(iPt) + "_ct" + std::to_string(iCt);

      for(int iScore=0, nScores=npLowerCutsForCorrBg.size(); iScore<nScores; ++iScore) {
        TH1* hCorrBg = fileIn->Get<TH1>((hPath + "hMass_bkgSum_NPgt" + to_string_with_precision(npLowerCutsForCorrBg.at(iScore), 2)).c_str());
        if(hCorrBg == nullptr) throw std::runtime_error("hCorrBg == nullptr for " + hPath);
        if(rebinFactor != 1) hCorrBg->Rebin(rebinFactor);
        if(smoothFactor != 0) hCorrBg->Smooth(smoothFactor);
        hCorrBg->SetLineColor(palette.at(iScore%palette.size()));
        hCorrBg->SetMarkerColor(palette.at(iScore%palette.size()));
        hCorrBg->SetLineWidth(3);
        ccCorrBg.cd();
        std::string drawOption = iScore == 0 ? "" : "same";
        if(smoothFactor != 0) drawOption += " HIST L";
        hCorrBg->Draw(drawOption.c_str());
        if(iCt==0) leg.AddEntry(hCorrBg, ("bdt(NP) > " + to_string_with_precision(npLowerCutsForCorrBg.at(iScore), 2)).c_str(), "L");
      } // npLowerCutsForCorrBg
      legTextCommon.Draw("same");
      leg.Draw("same");

      fileOut->cd();
      ccBoth.Write(("ccBoth" + ccNameSuffix).c_str());
      ccCorrBg.Write(("ccCorrBg" + ccNameSuffix).c_str());

      for(const auto& cc : {&ccCorrBg, &ccBoth}) {
        cc->Print((static_cast<std::string>(cc->GetName()) + ccNameSuffix + ".pdf").c_str(), "pdf");
      }
    } // ctRanges
  } // ptRanges
  fileOut->Close();
}
