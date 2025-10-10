namespace Particles {
enum Particles : short {
  kProton = 0,
  kPion,
  kElectorn,
  kKaon,
  kAll,
  nParticles
};
};

const std::array<std::string, Particles::nParticles> particleNames{"proton", "pion", "electorn", "kaon", "all"};

std::vector<std::string> GetDFNames(const std::string& fileName);
float calcBetaFromBetaGamma(float betaGamma);
short GetParicleByfPidIndex(UChar_t fPidIndex);

void tpc_qa(const std::string& fileName) {
  const std::vector<std::string> dirNames = GetDFNames(fileName);

  UChar_t fPidIndex;
  Float_t fTPCInnerParam;
  Float_t fBetaGamma;
  Float_t fTPCSignal;
  Float_t fNSigTPC;
  Float_t fNSigTOF;

  std::array<TH2D*, Particles::nParticles> hPdEdx;
  std::array<TH1D*, Particles::nParticles> hNSigmaTpc;
  std::array<TH2D*, Particles::nParticles> hPNSigmaTpc;
  std::array<TH1D*, Particles::nParticles> hNSigmaTof;
  std::array<TH2D*, Particles::nParticles> hPNSigmaTof;


  const int nBinsP = 100;
  const double lowP = 0.1;
  const double hiP = 10;

  std::vector<double> binEdgesP(nBinsP+1);
  const double logLowP = std::log10(lowP);
  const double logHiP = std::log10(hiP);
  const double logStepP = (logHiP - logLowP) / nBinsP;
  for(int iBin=0; iBin<=nBinsP; ++iBin) {
    binEdgesP.at(iBin) = std::pow(10., logLowP + iBin*logStepP);
  }

  const int nBinsDedx = 100;
  const double lowDedx = 0;
  const double hiDedx = 150;

  const int nBinsNsigma = 100;
  const double lowNsigma = -10;
  const double hiNsigma = 10;

  for(int kParticle=0; kParticle<Particles::nParticles; ++kParticle) {
    hPdEdx.at(kParticle) = new TH2D(("hPdEdx_" + particleNames.at(kParticle)).c_str(), particleNames.at(kParticle).c_str(), nBinsP, binEdgesP.data(), nBinsDedx, lowDedx, hiDedx);
    hPdEdx.at(kParticle)->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
    hPdEdx.at(kParticle)->GetYaxis()->SetTitle("dE/dx (a.u.)");
    hPdEdx.at(kParticle)->GetZaxis()->SetTitle("Entries");

    hNSigmaTpc.at(kParticle) = new TH1D(("hNSigmaTpc_" + particleNames.at(kParticle)).c_str(), particleNames.at(kParticle).c_str(), nBinsNsigma, lowNsigma, hiNsigma);
    hNSigmaTpc.at(kParticle)->GetXaxis()->SetTitle("N#sigma TPC");
    hNSigmaTpc.at(kParticle)->GetYaxis()->SetTitle("Entries");

    hPNSigmaTpc.at(kParticle) = new TH2D(("hPNSigmaTpc_" + particleNames.at(kParticle)).c_str(), particleNames.at(kParticle).c_str(), nBinsP, binEdgesP.data(), nBinsNsigma, lowNsigma, hiNsigma);
    hPNSigmaTpc.at(kParticle)->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
    hPNSigmaTpc.at(kParticle)->GetYaxis()->SetTitle("N#sigma TPC");
    hPNSigmaTpc.at(kParticle)->GetZaxis()->SetTitle("Entries");

    hNSigmaTof.at(kParticle) = new TH1D(("hNSigmaTof_" + particleNames.at(kParticle)).c_str(), particleNames.at(kParticle).c_str(), nBinsNsigma, lowNsigma, hiNsigma);
    hNSigmaTof.at(kParticle)->GetXaxis()->SetTitle("N#sigma TOF");
    hNSigmaTof.at(kParticle)->GetYaxis()->SetTitle("Entries");

    hPNSigmaTof.at(kParticle) = new TH2D(("hPNSigmaTof_" + particleNames.at(kParticle)).c_str(), particleNames.at(kParticle).c_str(), nBinsP, binEdgesP.data(), nBinsNsigma, lowNsigma, hiNsigma);
    hPNSigmaTof.at(kParticle)->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
    hPNSigmaTof.at(kParticle)->GetYaxis()->SetTitle("N#sigma TOF");
    hPNSigmaTof.at(kParticle)->GetZaxis()->SetTitle("Entries");
  }

  for(const auto& dirName : dirNames) {
    TFile* fileIn = TFile::Open(fileName.c_str());
    TTree* treeIn = fileIn->Get<TTree>((dirName + "/O2tpcskimv0tree").c_str());
    if(treeIn == nullptr) treeIn = fileIn->Get<TTree>((dirName + "/O2tpcskimv0wde").c_str());
    treeIn->SetBranchAddress("fPidIndex", &fPidIndex);
    treeIn->SetBranchAddress("fTPCInnerParam", &fTPCInnerParam);
    treeIn->SetBranchAddress("fBetaGamma", &fBetaGamma);
    treeIn->SetBranchAddress("fTPCSignal", &fTPCSignal);
    treeIn->SetBranchAddress("fNSigTPC", &fNSigTPC);
    treeIn->SetBranchAddress("fNSigTOF", &fNSigTOF);

    const int nEntries = treeIn->GetEntries();
    for(int iEntry=0; iEntry<nEntries; ++iEntry) {
      treeIn->GetEntry(iEntry);
      const float p = fTPCInnerParam;
      const float beta = calcBetaFromBetaGamma(fBetaGamma);
      const float dEdx = fTPCSignal;
      const float nSigmaTpc = fNSigTPC;
      const float nSigmaTof = fNSigTOF;

      const short particleId = GetParicleByfPidIndex(fPidIndex);

      hPdEdx.at(particleId)->Fill(p, dEdx);
      hPdEdx.at(Particles::kAll)->Fill(p, dEdx);

      hNSigmaTpc.at(particleId)->Fill(nSigmaTpc);
      hNSigmaTpc.at(Particles::kAll)->Fill(nSigmaTpc);

      hPNSigmaTpc.at(particleId)->Fill(p, nSigmaTpc);
      hPNSigmaTpc.at(Particles::kAll)->Fill(p, nSigmaTpc);

      hNSigmaTof.at(particleId)->Fill(nSigmaTof);
      hNSigmaTof.at(Particles::kAll)->Fill(nSigmaTof);

      hPNSigmaTof.at(particleId)->Fill(p, nSigmaTof);
      hPNSigmaTof.at(Particles::kAll)->Fill(p, nSigmaTof);
    }

    fileIn->Close();
  }

  TFile* fileOut = TFile::Open("tpc_qa.root", "recreate");
  for(int kParticle=0; kParticle<Particles::nParticles; ++kParticle) {
    hPdEdx.at(kParticle)->Write();
    hNSigmaTpc.at(kParticle)->Write();
    hPNSigmaTpc.at(kParticle)->Write();
    hNSigmaTof.at(kParticle)->Write();
    hPNSigmaTof.at(kParticle)->Write();
  }
  fileOut->Close();
}

std::vector<std::string> GetDFNames(const std::string& fileName) {
  TFile* fileIn = TFile::Open(fileName.c_str());
  if(fileIn == nullptr) {
    throw std::runtime_error("fileIn == nullptr");
  }

  std::vector<std::string> result;
  auto lok = fileIn->GetListOfKeys();
  for(const auto& k : *lok) {
    const std::string dirname = k->GetName();
    if(dirname == "parentFiles") continue;
    result.emplace_back(dirname);
  }
  fileIn->Close();

  return result;
}

float calcBetaFromBetaGamma(float betaGamma) {
  return betaGamma / std::sqrt(1 + betaGamma*betaGamma);
}

short GetParicleByfPidIndex(UChar_t fPidIndex) {
  switch(fPidIndex) {
    case static_cast<UChar_t>(0): return Particles::kElectorn;
    case static_cast<UChar_t>(2): return Particles::kPion;
    case static_cast<UChar_t>(3): return Particles::kKaon;
    case static_cast<UChar_t>(4): return Particles::kProton;
    default: return -1;
  }
}
