int SB_Statys2Type(int status);
void CD(TFile* file, std::string dir);

template <typename T>
void Vector2DResize(std::vector<std::vector<T>>& vec, int size_outer, int size_inner);

template <typename T>
void VectorResize(std::vector<T>& vec, int size);

enum eType : short {
  kPrompt = 0,
  kNonPrompt,
  kImpossible,
  kNumberOfTypes
};

void res_and_pull_qa(const std::string& fileName, int selectionFlag=1) {
  TFile* fileIn = TFile::Open(fileName.c_str());
  if(fileIn == nullptr) {
    throw std::runtime_error("fileIn == nullptr");
  }

  TFile* fileOut = TFile::Open("res_and_pull_qa.root", "recreate");

  std::vector<std::string> dirTypes{"prompt", "nonprompt", "impossible"};
  std::vector<std::string> histTypes{"residual", "corr", "pull"};

  TH1D* hSBType = new TH1D("hSBType", "hSBType", 7, -2.5, 4.5);
  const char* sb_titles[7] = {"", "Impossible", "Background", "Prompt", "Non-prompt", "Wrong swap", ""};
  for (int i=1;i<=7;i++) {
    hSBType->GetXaxis()->SetBinLabel(i, sb_titles[i-1]);
  }

  struct Variable {
    std::string name_;
    std::string name_in_tree_rec_;
    std::string name_in_tree_mc_;
    std::string name_in_tree_error_;
    std::string title_;
    std::string unit_;
    int nbins_plain_;
    float xlow_plain_;
    float xup_plain_;
    int nbins_res_;
    float xlow_res_;
    float xup_res_;
    int nbins_pull_;
    float xlow_pull_;
    float xup_pull_;
  };

  const int nbins = 400;

  std::vector<Variable> vars {
    {"P",  /*KF*/"fP",  /*MC*/"fP",       /*KF*/"fDeltaP",  "p",     "GeV/c", nbins, 0,    12,  nbins, -2, 2,     nbins, -5,  5},
    {"Pt", /*KF*/"fPt", /*MC*/"fPt",      /*KF*/"fDeltaPt", "p_{T}", "GeV/c", nbins, 0,    12,  nbins, -2, 2,     nbins, -5,  5},
    {"X",  /*KF*/"fX",  /*MC*/"fDecayX",  /*KF*/"fErrX",    "X",     "cm",    nbins, -0.5, 0.5, nbins, -0.1, 0.1, nbins, -5,  5},
    {"Y",  /*KF*/"fY",  /*MC*/"fDecayY",  /*KF*/"fErrY",    "Y",     "cm",    nbins, -0.5, 0.5, nbins, -0.1, 0.1, nbins, -5,  5},
    {"Z",  /*KF*/"fZ",  /*MC*/"fDecayZ",  /*KF*/"fErrZ",    "Z",     "cm",    nbins, -0.5, 0.5, nbins, -0.1, 0.1, nbins, -5,  5},
    {"T",  /*KF*/"fT",  /*MC*/"fDecayT",  /*KF*/"fDeltaT",  "T",     "ps",    nbins, -1,   5,   nbins, -3, 3,     nbins, -10, 10},
  };

  std::vector<std::vector<TH1D*>> hres;
  std::vector<std::vector<TH2D*>> hcorr;
  std::vector<std::vector<TH1D*>> hpull;
  Vector2DResize(hres, vars.size(), kNumberOfTypes);
  Vector2DResize(hcorr, vars.size(), kNumberOfTypes);
  Vector2DResize(hpull, vars.size(), kNumberOfTypes);

  for(int iVar=0; iVar<vars.size(); iVar++) {
    for(int iSBType=0; iSBType<kNumberOfTypes; iSBType++) {
      hres.at(iVar).at(iSBType) = new TH1D(("hRes_" + vars.at(iVar).name_ + "_" + dirTypes.at(iSBType)).c_str(), ("hRes_" + vars.at(iVar).name_).c_str(), vars.at(iVar).nbins_res_, vars.at(iVar).xlow_res_, vars.at(iVar).xup_res_);
      hcorr.at(iVar).at(iSBType) = new TH2D(("hCorr_" + vars.at(iVar).name_ + "_" + dirTypes.at(iSBType)).c_str(), ("hCorr_" + vars.at(iVar).name_).c_str(), vars.at(iVar).nbins_plain_, vars.at(iVar).xlow_plain_, vars.at(iVar).xup_plain_, vars.at(iVar).nbins_plain_, vars.at(iVar).xlow_plain_, vars.at(iVar).xup_plain_);
      hpull.at(iVar).at(iSBType) = new TH1D(("hPull_" + vars.at(iVar).name_ + "_" + dirTypes.at(iSBType)).c_str(), ("hPull_" + vars.at(iVar).name_).c_str(), vars.at(iVar).nbins_pull_, vars.at(iVar).xlow_pull_, vars.at(iVar).xup_pull_);

      hres.at(iVar).at(iSBType)->GetXaxis()->SetTitle((vars.at(iVar).title_ + "^{rec} - " + vars.at(iVar).title_ + "^{mc}, " + vars.at(iVar).unit_).c_str());
      hres.at(iVar).at(iSBType)->GetYaxis()->SetTitle("Entries");

      hcorr.at(iVar).at(iSBType)->GetXaxis()->SetTitle((vars.at(iVar).title_ + "^{mc}, " + vars.at(iVar).unit_).c_str());
      hcorr.at(iVar).at(iSBType)->GetYaxis()->SetTitle((vars.at(iVar).title_ + "^{rec}, " + vars.at(iVar).unit_).c_str());

      hpull.at(iVar).at(iSBType)->GetXaxis()->SetTitle(("(" + vars.at(iVar).title_ + "^{rec} - " + vars.at(iVar).title_ + "^{mc}) / #sigma_{" + vars.at(iVar).title_ + "^{rec}}").c_str());
      hpull.at(iVar).at(iSBType)->GetYaxis()->SetTitle("Entries");
    }
  }

  std::vector<float> value_mc, value_rec, value_err;
  VectorResize(value_mc, vars.size());
  VectorResize(value_rec, vars.size());
  VectorResize(value_err, vars.size());
  int sb_status, is_selected;

  auto lok = fileIn->GetListOfKeys();
  fileIn->cd();
  for(const auto& k : *lok) {
    const std::string dirname = k->GetName();
    if(dirname == "parentFiles") continue;

    TTree* treeKF = fileIn->Get<TTree>((dirname + "/O2hfcandlckf").c_str());
    TTree* treeMC = fileIn->Get<TTree>((dirname + "/O2hfcandlcmc").c_str());
    treeKF->SetBranchAddress("fSigBgStatus", &sb_status);
    treeKF->SetBranchAddress("fIsSelected", &is_selected);
    for(int iVar=0; iVar<vars.size(); iVar++) {
      treeKF->SetBranchAddress(vars.at(iVar).name_in_tree_rec_.c_str(), &value_rec.at(iVar));
      treeKF->SetBranchAddress(vars.at(iVar).name_in_tree_error_.c_str(), &value_err.at(iVar));
      treeMC->SetBranchAddress(vars.at(iVar).name_in_tree_mc_.c_str(), &value_mc.at(iVar));
    }

    const int Nentries = treeKF->GetEntries();
    for(int iEntry=0; iEntry<Nentries; iEntry++) {
      treeKF->GetEntry(iEntry);
      treeMC->GetEntry(iEntry);

      if(is_selected<selectionFlag) continue;

      hSBType->Fill(sb_status);
      const int sb_histotype = SB_Statys2Type(sb_status);

      for(int iVar=0; iVar<vars.size(); iVar++) {
        hres.at(iVar).at(sb_histotype)->Fill(value_rec.at(iVar) - value_mc.at(iVar));
        hcorr.at(iVar).at(sb_histotype)->Fill(value_mc.at(iVar), value_rec.at(iVar));
        hpull.at(iVar).at(sb_histotype)->Fill((value_rec.at(iVar) - value_mc.at(iVar)) / value_err.at(iVar));
      }
    } // iEntry
  } // GetListOfKeys

  fileOut->cd();
  hSBType->Write();
  for(int kType=0; kType<kNumberOfTypes; kType++) {
    for(int iVar=0; iVar<vars.size(); iVar++) {
      CD(fileOut, (histTypes.at(0) + "/" + dirTypes.at(kType)).c_str());
      hres.at(iVar).at(kType)->Write();
      CD(fileOut, (histTypes.at(1) + "/" + dirTypes.at(kType)).c_str());
      hcorr.at(iVar).at(kType)->Write();
      CD(fileOut, (histTypes.at(2) + "/" + dirTypes.at(kType)).c_str());
      hpull.at(iVar).at(kType)->Write();
    }
  }

  fileOut->Close();
  fileIn->Close();
}


int SB_Statys2Type(int status) {
  if      (status == 0) return kImpossible;
  else if (status == 1) return kPrompt;
  else if (status == 2) return kNonPrompt;
  else if (status == 3) return kImpossible;
  else return kImpossible;
}

void CD(TFile* file, std::string dir) {
  if(file == nullptr) throw std::runtime_error("CD() - file is nullptr");

  if(file->GetDirectory(dir.c_str()) == nullptr) file->mkdir(dir.c_str());
  file->cd(dir.c_str());
}

template <typename T>
void Vector2DResize(std::vector<std::vector<T>>& vec, int size_outer, int size_inner) {
  vec.resize(size_outer);
  for(auto& v : vec) {
    v.resize(size_inner);
  }
}

template <typename T>
void VectorResize(std::vector<T>& vec, int size) {
  vec.resize(size);
}
