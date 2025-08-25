const std::vector<std::string> particleNames{"proton", "pion", "electorn", "kaon", "all"};

void tpc_qa2(const std::string& fileName) {
  TFile* fileIn = TFile::Open(fileName.c_str(), "read");
  for(const auto& paName : particleNames) {
    const std::string priBra = particleNames.size() == 1 ? "" : paName == particleNames.front() ? "(" : paName == particleNames.back() ? ")" : "";
    TH1* h = fileIn->Get<TH1>(("hPdEdx_" + paName).c_str());
    TCanvas cc("cc", "");
    cc.SetCanvasSize(1200, 800);
    h->Draw("colz");
    cc.Print(("hPdEdx.pdf" + priBra).c_str(), "pdf");
  }
  fileIn->Close();
}
