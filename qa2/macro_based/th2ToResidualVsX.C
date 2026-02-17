int FindWhoseTH1LowerEdge(TAxis* histo, double value);

void th2ToResidualVsX(const std::string& fileName, const std::string& histoName) {
  gROOT->Macro("/home/oleksii/alidir/macros_on_git/qa2/exe_based/styles/mc_qa2.style.cc");

  std::vector<double> binEdges;
  for(int iStep=0; iStep<=40; ++iStep) {
    binEdges.emplace_back(iStep * 0.05);
  }

  TFile* fileIn = TFile::Open(fileName.c_str(), "read");
  if(fileIn == nullptr) throw std::runtime_error("File " + fileName + " is missing");

  TH2* histoIn = fileIn->Get<TH2>(histoName.c_str());
  if(histoIn == nullptr) throw std::runtime_error("Histogram " + histoName + " is missing");

  TGraphErrors grRes;
  grRes.SetName("grRes");
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

  grRes.SaveAs("grRes.root");

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
