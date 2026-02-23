int FindWhoseTH1LowerEdge(TAxis* histo, double value);

void th2ToResidualVsX(const std::string& fileName, const std::string& histoName) {
  gROOT->Macro("/home/oleksii/alidir/macros_on_git/qa2/exe_based/styles/mc_qa2.style.cc");

  std::vector<double> binEdges;
  double edge{0.};
  while(edge < 1.) {
    binEdges.emplace_back(edge);
    edge += 0.05;
  }
  while(edge <= 2.) {
    binEdges.emplace_back(edge);
    edge += 0.2;
  }

  TFile* fileIn = TFile::Open(fileName.c_str(), "read");
  if(fileIn == nullptr) throw std::runtime_error("File " + fileName + " is missing");

  TH2* histoIn = fileIn->Get<TH2>(histoName.c_str());
  if(histoIn == nullptr) throw std::runtime_error("Histogram " + histoName + " is missing");

  TGraphErrors grMean;
  TGraphErrors grSigma;
  grMean.GetXaxis()->SetTitle(histoIn->GetXaxis()->GetTitle());
  grMean.GetYaxis()->SetTitle(("#mu of "s + histoIn->GetYaxis()->GetTitle()).c_str());
  grSigma.GetXaxis()->SetTitle(histoIn->GetXaxis()->GetTitle());
  grSigma.GetYaxis()->SetTitle(("#sigma of "s + histoIn->GetYaxis()->GetTitle()).c_str());

  TAxis* xAxis = histoIn->GetXaxis();

  for(int iRange=0, nRanges=binEdges.size()-1; iRange<nRanges; ++iRange) {
    const double lo = binEdges.at(iRange);
    const double hi = binEdges.at(iRange+1);
    const int binLo = FindWhoseTH1LowerEdge(xAxis, lo);
    const int binHi = FindWhoseTH1LowerEdge(xAxis, hi) - 1;
    auto* histoProj = histoIn->ProjectionY("projY", binLo, binHi);
    const double grX = (binEdges.at(iRange+1) + binEdges.at(iRange)) / 2;
    const double grMeanY = histoProj->GetMean();
    const double grMeanYErr = histoProj->GetMeanError();
    const double grSigmaY = histoProj->GetStdDev();
    const double grSigmaYErr = histoProj->GetStdDevError();
    grMean.AddPoint(grX, grMeanY);
    grMean.SetPointError(grMean.GetN()-1, 0, grMeanYErr);
    grSigma.AddPoint(grX, grSigmaY);
    grSigma.SetPointError(grSigma.GetN()-1, 0, grSigmaYErr);
  }

  TFile* fileOut = TFile::Open("grRes.root", "recreate");
  grMean.Write("mean");
  grSigma.Write("sigma");
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
