void pt_NPweight_builder(const std::string& fileGenName, const std::string& fileFonllName, const std::string& fileInOutName) {
  TFile* fileGen = TFile::Open(fileGenName.c_str(), "read");
  TFile* fileFonll = TFile::Open(fileFonllName.c_str(), "read");

  TH1* histoGen = fileGen->Get<TH1>("hPtLambdaB");
  TH1* histoFonll = fileFonll->Get<TH1>("hPtFONLLBcent");

  histoGen->Rebin(10);

  const double integralGen = histoGen->Integral();
  const double integralFonll = histoFonll->Integral();

  histoGen->Scale(1./integralGen);
  histoFonll->Scale(1./integralFonll);

  const bool ok = histoFonll->Divide(histoGen);
  if(!ok) throw std::runtime_error("histoFonll->Divide(histoGen) was not ok");

  TFile* fileOut = TFile::Open(fileInOutName.c_str(), "update");
  histoFonll->Write("histoNPWeight");
  fileOut->Close();

  fileFonll->Close();
  fileGen->Close();
}
