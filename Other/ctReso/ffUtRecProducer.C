double getChi2BetweenTwoHistos(const TH1* hObs, const TH1* hRef);
void zeroNonDiagonalBins(TH2* histo);

void ffUtRecProducer(int seed = -1) {
  if(seed >= 0) gRandom->SetSeed(seed);

  const std::string fileName = "/home/oleksii/alidir/working/RMatrix/HF_LHC24h1b_All.757856/effs/efficiency_summary.root";
  const std::string histoName = "effs/prompt/pT_3_20/r_NPgt0.20_W";
  const double tau{0.206};
  const int nFills{1000000};
//   const std::vector<double> edges{0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8};

  TFile* fileResponse = TFile::Open(fileName.c_str());
  TH2* histoResponse = fileResponse->Get<TH2>(histoName.c_str());
  if(histoResponse == nullptr) throw std::runtime_error("histoResponse == nullptr");
//   zeroNonDiagonalBins(histoResponse);

  const int nBins = histoResponse->GetNbinsX();

  TH1* hEffSim = histoResponse->ProjectionY();

  std::vector<TH1*> histoMigration(nBins+1, nullptr);
  for(int iBin=1; iBin<=nBins; ++iBin) {
    histoMigration.at(iBin) = histoResponse->ProjectionX(("projX_" + std::to_string(iBin)).c_str(), iBin, iBin);
  }

  TH1* hExpo = dynamic_cast<TH1*>(histoResponse->ProjectionX()->Clone());
  hExpo->SetDirectory(nullptr);
  hExpo->Reset();
  hExpo->Sumw2(false);
  for(int iBin=1, nBins=hExpo->GetNbinsX(); iBin<=nBins; ++iBin) {
    const double lo = hExpo->GetBinLowEdge(iBin);
    const double hi = hExpo->GetBinLowEdge(iBin + 1);
    const double value = nFills * (std::exp(-lo/tau) - std::exp(-hi/tau));
    hExpo->SetBinContent(iBin, value);
  }

  TH1* hGen = dynamic_cast<TH1*>(histoResponse->ProjectionX()->Clone());
  hGen->SetDirectory(nullptr);
  hGen->Reset();
  hGen->Sumw2(false);

  TH1* hRecSampled = dynamic_cast<TH1*>(histoResponse->ProjectionX()->Clone());
  hRecSampled->SetDirectory(nullptr);
  hRecSampled->Reset();
  hRecSampled->Sumw2(false);

  TH1* hRecMultiplied = dynamic_cast<TH1*>(histoResponse->ProjectionX()->Clone());
  hRecMultiplied->SetDirectory(nullptr);
  hRecMultiplied->Reset();
  hRecMultiplied->Sumw2(false);

  for(int iFill=0; iFill<nFills; ++iFill) {
//     if(iFill%(nFills/20) == 0) std::cout << "iFill = " << iFill << "\n";
    const double tGen = gRandom->Exp(tau);
    hGen->Fill(tGen);
    const int binGen = hEffSim->FindBin(tGen);
    if(binGen<1 || binGen>nBins) continue;
    if(gRandom->Uniform(1) > hEffSim->GetBinContent(binGen)) continue;
    hRecSampled->Fill(histoMigration.at(binGen)->GetRandom());
  }
//   hRecSampled = dynamic_cast<TH1*>(hRecSampled->Rebin(edges.size() - 1, hRecSampled->GetName(), edges.data()));

  for(int iBinRec=1; iBinRec<=nBins; ++iBinRec) {
    double result{};
    for(int iBinGen=1; iBinGen<=nBins; ++iBinGen) {
      result += histoResponse->GetBinContent(iBinRec, iBinGen) * hGen->GetBinContent(iBinGen);
    }
    hRecMultiplied->SetBinContent(iBinRec, result);
  }

  hExpo->SetName("hExpo");
  hGen->SetName("hGen");
  hRecSampled->SetName("hCorrYieldsPrompt");
  hRecMultiplied->SetName("hRecMultiplied");

  std::cout << "hGen vs hExpo chi2 / ndf = " << getChi2BetweenTwoHistos(hGen, hExpo) << " / " << hGen->GetNbinsX() << "\t\t";
  std::cout << "hRecSampled vs hRecMultiplied chi2 / ndf = " << getChi2BetweenTwoHistos(hRecSampled, hRecMultiplied) << " / " << hRecSampled->GetNbinsX() << "\n";


  TFile* fileOut = TFile::Open(("unitTestRec.nFills" + std::to_string(nFills) + ".seed" + std::to_string(seed) + ".root").c_str(), "recreate");
  hExpo->Write();
  hGen->Write();
  hRecSampled->Write();
  hRecMultiplied->Write();
  fileOut->Close();
}

double getChi2BetweenTwoHistos(const TH1* hObs, const TH1* hRef) {
  double chi2{};
  for(int iBin=1, nBins=hRef->GetNbinsX(); iBin<=nBins; ++iBin) {
    const double diff = hObs->GetBinContent(iBin) - hRef->GetBinContent(iBin);
    const double err = hObs->GetBinError(iBin);
    chi2 += diff*diff / err / err;
  }

  return chi2;
}

void zeroNonDiagonalBins(TH2* histo) {
  const int nBins = histo->GetNbinsX();
  for(int iBin=1; iBin<=nBins; ++iBin) {
    for(int jBin=1; jBin<=nBins; ++jBin) {
      if(iBin != jBin) histo->SetBinContent(iBin, jBin, 0.);
    }
  }
}
