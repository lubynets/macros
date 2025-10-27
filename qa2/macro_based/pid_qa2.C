const std::vector<std::string> PidDetectors{"Tpc", "Tof", "TpcTof"};
const std::vector<std::string> ProngSpecies{"Pr", "Ka", "Pi"};

void pid_qa2(const std::string& fileName) {
  TFile* fileIn = TFile::Open(fileName.c_str(), "read");
  for(int iPs=0, nPss=ProngSpecies.size(); iPs<nPss; iPs++) {
    for(int iDet=0, nDets=PidDetectors.size(); iDet<nDets; iDet++) {
      TH2* h = fileIn->Get<TH2>(("PidQa/h2nSig" + PidDetectors.at(iDet) + ProngSpecies.at(iPs) + "VsPt").c_str());
      TCanvas cc("cc", "");
      cc.SetCanvasSize(1200, 800);
      cc.SetRightMargin(0.12);
      cc.SetLogx();
      cc.SetLogz();
      h->GetZaxis()->CenterTitle();
      h->Draw("colz");
      const std::string priBra = iPs+iDet == 0 ? "(" : iPs*iDet == (nPss-1)*(nDets-1) ? ")" : "";
      cc.Print(("pid_qa.pdf" + priBra).c_str(), "pdf");
    } // nPidDetectors
  } // nProngSpecies

  fileIn->Close();
}
