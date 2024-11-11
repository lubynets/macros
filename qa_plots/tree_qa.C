void tree_qa(const std::string& fileName = "/home/oleksii/alidir/sandbox/AnalysisResults_trees.nolite.root", const std::string& treeName="O2hfcand3pbase") {

  TFile* fileIn = TFile::Open(fileName.c_str());
  if(fileIn == nullptr) {
    throw std::runtime_error("fileIn == nullptr");
  }

  TFile* fileOut = TFile::Open("fileOut.root", "recreate");

  float posx;
  TH1F* hposx = new TH1F("hposx", "HPOSX", 100, -0.1, 0.1);

  auto lok = fileIn->GetListOfKeys();
  fileIn->cd();
  for(const auto& k : *lok) {
    const std::string dirname = k->GetName();
    if(dirname == "parentFiles") continue;
    TTree* treeIn = fileIn->Get<TTree>((dirname + "/" + treeName).c_str());

    treeIn->SetBranchAddress("fPosX", &posx);

    const int Nentries = treeIn->GetEntries();
    for(int iEntry=0; iEntry<Nentries; iEntry++) {
      treeIn->GetEntry(iEntry);
      hposx->Fill(posx);
    }
  }


  fileOut->cd();

  hposx->Write();

  fileOut->Close();
  fileIn->Close();
}
