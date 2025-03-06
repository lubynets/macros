std::vector<std::string> GetDFNames(const std::string& fileName);

void IdentEntriesChecker(const std::string& fileName, int nDF=1) {
  auto dirNames = GetDFNames(fileName);

  const int nbins = 100000;

  TH1D* hP = new TH1D("hP", "", nbins, 0, 10);
  TH1D* hPt = new TH1D("hPt", "", nbins, 0, 10);
  TH1D* hY = new TH1D("hY", "", nbins, -1, 1);
  TH1D* hMPKPi = new TH1D("hMPKPi", "", nbins, 1.7, 2.7);

  int iDF{0};
  for(auto& dir : dirNames) {
    if(iDF > 0) break;
    std::cout << "dirname = " << dir << "\n";
    TFile* fileIn = TFile::Open(fileName.c_str(), "read");
    TTree* treeKF = fileIn->Get<TTree>((dir + "/O2hfcandlckf").c_str());
    TTree* treeLite = fileIn->Get<TTree>((dir + "/O2hfcandlclite").c_str());

    float P, Pt, Y, M;
    treeKF->SetBranchAddress("fP", &P);
    treeKF->SetBranchAddress("fPt", &Pt);
    treeLite->SetBranchAddress("fY", &Y);
    treeKF->SetBranchAddress("fMassInv", &M);

    if(treeKF->GetEntries() != treeLite->GetEntries()) throw std::runtime_error("treeKF->GetEntries() != treeLite->GetEntries()");

    for(int iEntry=0; iEntry<treeKF->GetEntries(); iEntry++) {
      treeKF->GetEntry(iEntry);
      treeLite->GetEntry(iEntry);

      hP->Fill(P);
      hPt->Fill(Pt);
      hY->Fill(Y);
      hMPKPi->Fill(M);

    }
    fileIn->Close();
    iDF++;
  }
  TFile* fileOut = TFile::Open("fileOut.root", "recreate");
  hP->Write();
  hPt->Write();
  hY->Write();
  hMPKPi->Write();
  fileOut->Close();
}

std::vector<std::string> GetDFNames(const std::string& fileName) {
  TFile* fileIn = TFile::Open(fileName.c_str(), "read");
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
