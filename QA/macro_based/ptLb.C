std::vector<std::string> GetDFNames(const std::string& fileName);

void pt_gen_builder_alitree(const std::string& fileListName) {
  TH1D* hLbPt = new TH1D("hLbPt", "hLbPt", 2000, 0, 20);
  hLbPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hLbPt->GetYaxis()->SetTitle("Entries");

  std::ifstream fileList(fileListName);

    if (!fileList) {
        std::cerr << "Error: Could not open file \"" << fileListName << "\" for reading.\n";
        return;
    }

    std::string fileName;
    while (std::getline(fileList, fileName)) {

      std::cout << "Processing file " << fileName << "\n";
      const std::vector<std::string> dirNames = GetDFNames(fileName);

      for(const auto& dirName : dirNames) {
        TFile* fileIn = TFile::Open(fileName.c_str());
        TTree* treeIn = fileIn->Get<TTree>((dirName + "/O2mcparticle_001").c_str());
        std::cout << "Processing DF " << dirName << "\n";

        TH1D* h1 = new TH1D("h1", "h1", 2000, 0, 20);
        treeIn->Draw("TMath::Sqrt(fPx*fPx + fPy*fPy)>>h1", "fPdgCode == 5122");
        hLbPt->Add(h1);

        delete h1;
        fileIn->Close();
      }
    }

  TFile* fileOut = TFile::Open("fileOut.root", "recreate");
  hLbPt->Write();
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
    if(dirname.substr(0, 3) != "DF_") continue;
    result.emplace_back(dirname);
  }
  fileIn->Close();

  return result;
}
