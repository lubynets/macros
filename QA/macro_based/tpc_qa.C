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

  std::array<TH2D*, Particles::nParticles> hPdEdx;

  for(int kParticle=0; kParticle<Particles::nParticles; ++kParticle) {
    hPdEdx.at(kParticle) = new TH2D(("hPdEdx_" + particleNames.at(kParticle)).c_str(), particleNames.at(kParticle).c_str(), 100, 0.01, 10, 100, 0, 150);
    hPdEdx.at(kParticle)->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
    hPdEdx.at(kParticle)->GetYaxis()->SetTitle("dE/dx (a.u.)");
  }

  for(const auto& dirName : dirNames) {
    TFile* fileIn = TFile::Open(fileName.c_str());
    TTree* treeIn = fileIn->Get<TTree>((dirName + "/O2tpcskimv0tree").c_str());
    treeIn->SetBranchAddress("fPidIndex", &fPidIndex);
    treeIn->SetBranchAddress("fTPCInnerParam", &fTPCInnerParam);
    treeIn->SetBranchAddress("fBetaGamma", &fBetaGamma);
    treeIn->SetBranchAddress("fTPCSignal", &fTPCSignal);

    const int nEntries = treeIn->GetEntries();
    for(int iEntry=0; iEntry<nEntries; ++iEntry) {
      treeIn->GetEntry(iEntry);
      const float p = fTPCInnerParam;
      const float beta = calcBetaFromBetaGamma(fBetaGamma);
      const float dEdx = fTPCSignal;
      hPdEdx.at(GetParicleByfPidIndex(fPidIndex))->Fill(p, dEdx);
      hPdEdx.at(Particles::kAll)->Fill(p, dEdx);
    }

    fileIn->Close();
  }

  TFile* fileOut = TFile::Open("tpc_qa.root", "recreate");
  for(int kParticle=0; kParticle<Particles::nParticles; ++kParticle) {
    hPdEdx.at(kParticle)->Write();
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
