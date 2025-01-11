void CheckTreesNEntriesConsistensy(std::vector<TTree*> trees);
void GetEntryForTreesVector(std::vector<TTree*> trees, int iEntry);
void CD(TFile* file, std::string dir);

void mass_qa(const std::string& fileName, int selectionFlag=1) {

  const int LcFlag{2};
  const int PromptOrigin{1};
  const int NonPromptOrigin{2};

  TFile* fileIn = TFile::Open(fileName.c_str());
  if(fileIn == nullptr) {
    throw std::runtime_error("fileIn == nullptr");
  }

  TFile* fileOut = TFile::Open("mass_qa.root", "recreate");

  enum eType : short {
    kPrompt = 0,
    kNonPrompt,
    kBackground,
    kWrongSwap,
    kNumberOfTypes
  };

  std::vector<std::string> dirTypes{"prompt", "nonprompt", "background", "wrongswap"};

  std::vector<TH1F*> hmass;
  hmass.resize(kNumberOfTypes);
  int iType{0};
  for(auto& h : hmass) {
    h = new TH1F(("hMass_" + dirTypes.at(iType)).c_str(), "hMass", 600, 1.98, 2.58);
    h->GetXaxis()->SetTitle("m_{pK#pi}, GeV/c^{2}");
    h->GetYaxis()->SetTitle("Entries");
    ++iType;
  }

  float massPKPi, massPiKP;
  int isSelPKPi, isSelPiKP;
  Char_t flag_c, origin_c, swapped_c;

  auto lok = fileIn->GetListOfKeys();
  fileIn->cd();
  for(const auto& k : *lok) {
    const std::string dirname = k->GetName();
    if(dirname == "parentFiles") continue;

    TTree* treeBase = fileIn->Get<TTree>((dirname + "/O2hfcand3pbase").c_str());

    TTree* treeKF = fileIn->Get<TTree>((dirname + "/O2hfcand3pkf").c_str());
    treeKF->SetBranchAddress("fKfMassPKPi", &massPKPi);
    treeKF->SetBranchAddress("fKfMassPiKP", &massPiKP);

    TTree* treeSel = fileIn->Get<TTree>((dirname + "/O2hfsellc").c_str());
    treeSel->SetBranchAddress("fIsSelLcToPKPi", &isSelPKPi);
    treeSel->SetBranchAddress("fIsSelLcToPiKP", &isSelPiKP);

    TTree* treeMcRec = fileIn->Get<TTree>((dirname + "/O2hfcand3pmcrec").c_str());
    treeMcRec->SetBranchAddress("fFlagMcMatchRec", &flag_c);
    treeMcRec->SetBranchAddress("fOriginMcRec", &origin_c);
    treeMcRec->SetBranchAddress("fIsCandidateSwapped", &swapped_c);

    std::vector<TTree*> treesVec{treeBase, treeKF, treeSel, treeMcRec};

    CheckTreesNEntriesConsistensy(treesVec);
    const int Nentries = treeBase->GetEntries();
    for(int iEntry=0; iEntry<Nentries; iEntry++) {
      GetEntryForTreesVector(treesVec, iEntry);
      int flag = (int)flag_c;
      int origin = (int)origin_c;
      int swapped = (int)swapped_c;

      if(std::abs(flag) == LcFlag) { // is signal
        if(swapped == 0) { // PKPi is true signal
          if(isSelPKPi >= selectionFlag) {
            if(origin == PromptOrigin)         hmass.at(kPrompt)->Fill(massPKPi);
            else if(origin == NonPromptOrigin) hmass.at(kNonPrompt)->Fill(massPKPi);
          }
          if(isSelPiKP >= selectionFlag) hmass.at(kWrongSwap)->Fill(massPiKP);
        } else { // PiKP is true signal
          if(isSelPiKP >= selectionFlag) {
            if(origin == PromptOrigin)         hmass.at(kPrompt)->Fill(massPiKP);
            else if(origin == NonPromptOrigin) hmass.at(kNonPrompt)->Fill(massPiKP);
          }
          if(isSelPKPi >= selectionFlag) hmass.at(kWrongSwap)->Fill(massPKPi);
        }
      } else {
        if(isSelPKPi >= selectionFlag) hmass.at(kBackground)->Fill(massPKPi);
        if(isSelPiKP >= selectionFlag) hmass.at(kBackground)->Fill(massPiKP);
      }
    }
  }

  for(int kType=0; kType<kNumberOfTypes; kType++) {
    CD(fileOut, dirTypes.at(kType));
    hmass.at(kType)->Write();;
  }

  fileOut->Close();
  fileIn->Close();
}

void CheckTreesNEntriesConsistensy(std::vector<TTree*> trees) {
  std::vector<int> Nentries;
  for(auto& t : trees) {
    Nentries.push_back(t->GetEntries());
  }
  if(!(std::equal(Nentries.begin() + 1, Nentries.end(), Nentries.begin()))) {
    throw std::runtime_error("CheckTreesNEntriesConsistensy() failed!");
  }
}

void GetEntryForTreesVector(std::vector<TTree*> trees, int iEntry) {
  for(auto& t : trees) {
    t->GetEntry(iEntry);
  }
}

void CD(TFile* file, std::string dir) {
  if(file == nullptr) throw std::runtime_error("Helper::CD() - file is nullptr");

  if(file->GetDirectory(dir.c_str()) == nullptr) file->mkdir(dir.c_str());
  file->cd(dir.c_str());
}
