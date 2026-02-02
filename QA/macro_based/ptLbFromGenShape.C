int FindWhoseTH1LowerEdge(TAxis* histo, double value);

void ptLbFromGenShape(const std::string& fileName) {
  TFile* fileIn = TFile::Open(fileName.c_str(), "read");
  TH2* histoIn = fileIn->Get<TH2>("hf-task-mc-gen-pt-rap-shapes/BeautyHadrons/hRapVsPt5122");

  const double yMax{0.8};

  const int binYLow = FindWhoseTH1LowerEdge(histoIn->GetYaxis(), -yMax);
  const int binYUp = FindWhoseTH1LowerEdge(histoIn->GetYaxis(), yMax) - 1;

  std::cout << "binYLow = " << binYLow << ", binYUp = " << binYUp << "\n";

  TH1* histoOut = histoIn->ProjectionX("projX", binYLow, binYUp);

  TFile* fileOut = TFile::Open("fileOut.root", "recreate");
  histoOut->Write("hPtLambdaB");
  fileOut->Close();

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
