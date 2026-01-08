template<typename T>
inline std::string to_string_with_precision(const T a_value, const int n=2) {
  std::ostringstream out;
  out.precision(n);
  out << std::fixed << a_value;
  return out.str();
}

const std::vector<Color_t> palette{kBlue, kRed, kGreen+2, kMagenta, kBlack, kOrange+2};

void minvBuilder() {
  gStyle->SetOptStat(0);
  const std::string fileName = "/home/oleksii/alidir/working/cutVar/data/input/mass_bdt_qa_thn.HL.data.HF_LHC23_pass4_Thin_2P3PDstar.574294.ctbin1.BGwise.root";
  TFile* fileIn = TFile::Open(fileName.c_str(), "read");

  const std::vector<float> ptRanges = {0, 1, 2, 3, 4, 5, 8, 12, 20};
  const std::vector<float> ctRanges = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 5.0};
  const std::vector<std::vector<float>> bgScores{
    {0.01, 0.02, 0.03, 0.04, 0.05},
    {0.01, 0.02, 0.03, 0.04, 0.05},
    {0.01, 0.02, 0.03, 0.04, 0.05},
    {0.01, 0.02, 0.03, 0.04, 0.05},
    {0.01, 0.02, 0.03, 0.04, 0.05},
    {0.01, 0.02, 0.03, 0.04, 0.05},
    {0.02, 0.04, 0.06, 0.08, 0.10},
    {0.04, 0.08, 0.12, 0.16, 0.20}
  };

  for(int iPt=0, nPts=ptRanges.size()-1; iPt<nPts; ++iPt) {
    for(int iCt=0, nCts=ctRanges.size()-1; iCt<nCts; ++iCt) {
      std::vector<TH1*> histos(bgScores.at(iPt).size(), nullptr);

      for(int iScore=0, nScores=bgScores.at(iPt).size(); iScore<nScores; ++iScore) {
        const std::string hName = "pT_" + to_string_with_precision(ptRanges.at(iPt), 0) + "_" + to_string_with_precision(ptRanges.at(iPt+1), 0) + "/T_" + to_string_with_precision(ctRanges.at(iCt)) + "_" + to_string_with_precision(ctRanges.at(iCt+1)) + "/hM_BGlt" + to_string_with_precision(bgScores.at(iPt).at(iScore), 2);
        TH1* hMass = fileIn->Get<TH1>(hName.c_str());
        if(hMass == nullptr) throw std::runtime_error(hName.c_str());
        hMass->Rebin(4);
        hMass->SetLineColor(palette.at(iScore));
        hMass->SetLineWidth(2);
        hMass->SetMinimum(0);
        histos.at(iScore) = hMass;
      }

      TCanvas cc("cc", "");
      cc.SetCanvasSize(1000, 1000);
      for(int iScore=bgScores.at(iPt).size()-1; iScore>=0; --iScore) {
        const std::string drawOption = iScore == bgScores.at(iPt).size()-1 ? "PE" : "PE same";
        histos.at(iScore)->Draw(drawOption.c_str());
      }

      TPaveText legText(0.18, 0.78, 0.35, 0.88, "brNDC");
      legText.AddText(("#it{p}_{T} #in (" + to_string_with_precision(ptRanges.at(iPt), 0) + "; " + to_string_with_precision(ptRanges.at(iPt+1), 0) + ") GeV/#it{c}").c_str());
      legText.AddText(("t #in (" + to_string_with_precision(ctRanges.at(iCt)) + "; " + to_string_with_precision(ctRanges.at(iCt+1)) + ") ps").c_str());
      legText.SetFillColor(0);
      legText.SetTextSize(0.03);
      legText.SetTextFont(62);
      legText.Draw("same");

      cc.Print(("pt" + std::to_string(iPt) + "_ct" + std::to_string(iCt) + ".pdf").c_str(), "pdf");
    }
  }

  fileIn->Close();
}
