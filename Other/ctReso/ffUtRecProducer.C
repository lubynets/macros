void ffUtRecProducer() {
  const std::string fileName = "/home/oleksii/alidir/working/RMatrix/HF_LHC24h1b_All.757856/effs/efficiency_summary.root";
  const std::string histoName = "effs/prompt/pT_3_20/r_NPgt0.20_W";
  const double tau{0.206};
  const int nFills{10'000'000};
  const std::vector<double> edges{0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8};

  TFile* fileResponse = TFile::Open(fileName.c_str());
  TH2* histoResponse = fileResponse->Get<TH2>(histoName.c_str());
  if(histoResponse == nullptr) throw std::runtime_error("histoResponse == nullptr");

  const int nBins = histoResponse->GetNbinsX();

  TH1* hEffSim = histoResponse->ProjectionY();

  std::vector<TH1*> histoMigration(nBins+1, nullptr);
  for(int iBin=1; iBin<=nBins; ++iBin) {
    histoMigration.at(iBin) = histoResponse->ProjectionX(("projX_" + std::to_string(iBin)).c_str(), iBin, iBin);
  }

  TH1* hRec = dynamic_cast<TH1*>(histoResponse->ProjectionX()->Clone());
  hRec->SetDirectory(nullptr);
  hRec->Reset();
  hRec->Sumw2(false);

  for(int iFill=0; iFill<nFills; ++iFill) {
    if(iFill%(nFills/20) == 0) std::cout << "iFill = " << iFill << "\n";
    const double tGen = gRandom->Exp(tau);
    const int binGen = hEffSim->FindBin(tGen);
    if(binGen<1 || binGen>nBins) continue;
    if(gRandom->Uniform(1) > hEffSim->GetBinContent(binGen)) continue;
    hRec->Fill(histoMigration.at(binGen)->GetRandom());
  }
  hRec = dynamic_cast<TH1*>(hRec->Rebin(edges.size() - 1, hRec->GetName(), edges.data()));

  hRec->SetName("hCorrYieldsPrompt");

  hRec->SaveAs("unitTestRec.root");
}
