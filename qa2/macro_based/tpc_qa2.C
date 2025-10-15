const std::vector<std::string> particleNames{"proton", "pion", "electorn", "kaon", "all"};
const std::vector<double> particleMasses{0.938272, 0.139571, 0.000511, 0.493677};
const std::vector<std::string> histoNames{"hPdEdx", "hNSigmaTpc", "hPNSigmaTpc", "hNSigmaTof", "hPNSigmaTof", "hPNoMatchedTof"};

void tpc_qa2(const std::string& fileName, bool isParticleWise=true, bool isDrawBB=true) {
  gStyle->SetOptStat("nemruo");

  const std::array<double, 5> bbParams {
    0.17489300668239594,
    4.753570079803467,
    0.0036148501094430685,
    2.33624005317688,
    0.9984580278396606,
  };

  TF1* funcBB = new TF1("funcBB", "([0]*([1]- log([2]+pow(x*[5],-1.*[4])) - pow((x*[5]/sqrt(1+x*x*[5]*[5])),[3]))*50./pow((x*[5]/sqrt(1+x*x*[5]*[5])),[3]))*pow(1.0,2.3)", 0.1, 10);
  for(int iPar=0; iPar<5; ++iPar) {
    funcBB->SetParameter(iPar, bbParams.at(iPar));
  }

  std::vector<TF1*> funcBBParticle(particleMasses.size(), nullptr);
  for(int iParticle=0, nParticles=particleMasses.size(); iParticle<nParticles; ++iParticle) {
    funcBBParticle.at(iParticle) = dynamic_cast<TF1*>(funcBB->Clone());
    funcBBParticle.at(iParticle)->SetParameter(5, 1. / particleMasses.at(iParticle));
    funcBBParticle.at(iParticle)->SetLineWidth(3);
  }

  TFile* fileIn = TFile::Open(fileName.c_str(), "read");
  for(const auto& paName : particleNames) {
    for(const auto& hiName : histoNames) {
      std::string priBra;
      std::string ccName;
      if(isParticleWise) {
        priBra = histoNames.size() == 1 ? "" : hiName == histoNames.front() ? "(" : hiName == histoNames.back() ? ")" : "";
        ccName = paName;
      } else {
        priBra = particleNames.size() == 1 ? "" : paName == particleNames.front() ? "(" : paName == particleNames.back() ? ")" : "";
        ccName = hiName;
      }
      TH1* h = fileIn->Get<TH1>((hiName + "_" + paName).c_str());
      TCanvas cc("cc", "");
      cc.SetCanvasSize(1200, 800);
      cc.SetRightMargin(0.12);
      cc.SetLogx(std::string(h->GetXaxis()->GetTitle()) == "#it{p} (GeV/#it{c})");
      cc.SetLogz();
      h->GetZaxis()->CenterTitle();
      h->Draw("colz");
      if(isDrawBB && hiName == "hPdEdx") {
        for(int iParticle=0, nParticles=particleMasses.size(); iParticle<nParticles; ++iParticle) {
          funcBBParticle.at(iParticle)->Draw("same");
        }
      }
      cc.Print((ccName + ".pdf" + priBra).c_str(), "pdf");
    }
  }
  fileIn->Close();
}
