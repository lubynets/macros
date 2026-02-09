int FindWhoseTH1LowerEdge(TAxis* histo, double value);

void th2ToResidualVsX() {
  gROOT->Macro("/home/oleksii/alidir/macros_on_git/qa2/exe_based/styles/mc_qa2.style.cc");

  const std::string fileName{"/home/oleksii/alidir/working/ctReso/ct_mc_qa.HL.mc.HF_LHC24h1b_All.595984.root"};
  const std::string histoName{"timeDiffVsRec"};

  const std::array binEdges{0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8};

  TFile* fileIn = TFile::Open(fileName.c_str(), "read");
  TH2* histoIn = fileIn->Get<TH2>(histoName.c_str());

  TGraphErrors grRes;
  grRes.GetXaxis()->SetTitle(histoIn->GetXaxis()->GetTitle());
  grRes.GetYaxis()->SetTitle(("std. dev. of "s + histoIn->GetYaxis()->GetTitle()).c_str());

  TAxis* xAxis = histoIn->GetXaxis();

  for(int iRange=0, nRanges=binEdges.size()-1; iRange<nRanges; ++iRange) {
    const double lo = binEdges.at(iRange);
    const double hi = binEdges.at(iRange+1);
    const int binLo = FindWhoseTH1LowerEdge(xAxis, lo);
    const int binHi = FindWhoseTH1LowerEdge(xAxis, hi) - 1;
    auto* histoProj = histoIn->ProjectionY("projY", binLo, binHi);
    const double grX = (binEdges.at(iRange+1) + binEdges.at(iRange)) / 2;
    const double grY = histoProj->GetStdDev();
    const double grYErr = histoProj->GetStdDevError();
    grRes.AddPoint(grX, grY);
    grRes.SetPointError(grRes.GetN()-1, 0, grYErr);
  }

  TCanvas cc("cc", "", 1200, 800);
  grRes.Draw("APE");

  cc.Print("grRes.pdf", "pdf");

  fileIn->Close();
}

int FindWhoseTH1LowerEdge(TAxis* axis, double value) {
  const int approxBinNumber = axis->FindBin(value);
  for(int iBin = approxBinNumber-2; iBin<=approxBinNumber+2; ++iBin) {
    const double lowerEdge = axis->GetBinLowEdge(iBin);
    if(std::fabs(lowerEdge - value) < 1e-4) return iBin;
  }
  throw std::runtime_error("FindWhoseTH1LowerEdge() - no bin with lower edge " + std::to_string(value));
}
