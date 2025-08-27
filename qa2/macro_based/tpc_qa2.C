const std::vector<std::string> particleNames{"proton", "pion", "electorn", "kaon", "all"};
const std::vector<std::string> histoNames{"hPdEdx", "hNSigma", "hPNSigma"};

void tpc_qa2(const std::string& fileName, bool isParticleWise=true) {
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
      cc.SetLogx(std::string(h->GetXaxis()->GetTitle()) == "#it{p} (GeV/#it{c})");
      cc.SetLogz();
      h->Draw("colz");
      cc.Print((ccName + ".pdf" + priBra).c_str(), "pdf");
    }
  }
  fileIn->Close();
}
