std::vector<std::string> GetDFNames(const std::string& fileName);

void treeCounter(const std::string& fileName) {
  auto dirNames = GetDFNames(fileName);

  int nPromptKF{0};
  int nNonPromptKF{0};
  int nPromptGen{0};
  int nNonPromptGen{0};

  for(auto& dir : dirNames) {
    TFile* fileIn = TFile::Open(fileName.c_str(), "read");
    TTree* treeKF = fileIn->Get<TTree>((dir + "/O2hfcandlckf").c_str());
    TTree* treeLite = fileIn->Get<TTree>((dir + "/O2hfcandlclite").c_str());
    TTree* treeGen = fileIn->Get<TTree>((dir + "/O2hfcandlcfullp").c_str());
    int sbStatusKF;
    char sbStatusGen;
    float yKF, yGen;
    treeKF->SetBranchAddress("fSigBgStatus", &sbStatusKF);
    treeLite->SetBranchAddress("fY", &yKF);
    treeGen->SetBranchAddress("fOriginMcGen", &sbStatusGen);
    treeGen->SetBranchAddress("fY", &yGen);

    if(treeKF->GetEntries() != treeLite->GetEntries()) throw std::runtime_error("treeKF->GetEntries() != treeLite->GetEntries()");

    for(int iEntry=0; iEntry<treeKF->GetEntries(); iEntry++) {
      treeKF->GetEntry(iEntry);
      treeLite->GetEntry(iEntry);
      if(std::abs(yKF) > 0.8) continue;

      if(sbStatusKF == 1) ++nPromptKF;
      if(sbStatusKF == 2) ++nNonPromptKF;
    }

    for(int iEntry=0; iEntry<treeGen->GetEntries(); iEntry++) {
      treeGen->GetEntry(iEntry);
      if(std::abs(yGen) > 0.8) continue;

      if(sbStatusGen == 1) ++nPromptGen;
      if(sbStatusGen == 2) ++nNonPromptGen;
    }
    fileIn->Close();
  }

  std::cout << "nPromptKF = " << nPromptKF << "\n";
  std::cout << "nNonPromptKF = " << nNonPromptKF << "\n";
  std::cout << "nPromptGen = " << nPromptGen << "\n";
  std::cout << "nNonPromptGen = " << nNonPromptGen << "\n";
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
