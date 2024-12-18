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

void mc_qa(const std::string& fileName, int selectionFlag=0) {
  TFile* fileIn = TFile::Open(fileName.c_str());
  if(fileIn == nullptr) {
    throw std::runtime_error("fileIn == nullptr");
  }

  TFile* fileOut = TFile::Open("mc_qa.root", "recreate");

  std::vector<std::string> dirTypes{"prompt", "nonprompt", "impossible"};
  std::vector<std::string> histTypes{"mc", "rec", "residual", "corr", "pull"};

  TH1D* hSBType = new TH1D("hSBType", "hSBType", 7, -2.5, 4.5);
  const char* sb_titles[7] = {"", "Impossible", "Background", "Prompt", "Non-prompt", "Wrong swap", " "};
  for (int i=1;i<=7;i++) {
    hSBType->GetXaxis()->SetBinLabel(i, sb_titles[i-1]);
  }

  struct Variable {
    std::string name_;
    std::string name_in_tree_rec_;
    bool rec_from_kf_;
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
    {"P",   /*KF*/"fP",    true,  /*MC*/"fP",       /*KF*/"fDeltaP",  "p",      "GeV/c", nbins, 0,    12,  nbins, -1, 1,       nbins, -5,  5 },
    {"Pt",  /*KF*/"fPt",   true,  /*MC*/"fPt",      /*KF*/"fDeltaPt", "p_{T}",  "GeV/c", nbins, 0,    12,  nbins, -1, 1,       nbins, -5,  5 },
    {"Xsv", /*KF*/"fX",    true,  /*MC*/"fDecayX",  /*KF*/"fErrX",    "X_{SV}", "cm",    nbins, -0.3, 0.3, nbins, -0.05, 0.05, nbins, -5,  5 },
    {"Ysv", /*KF*/"fY",    true,  /*MC*/"fDecayY",  /*KF*/"fErrY",    "Y_{SV}", "cm",    nbins, -0.3, 0.3, nbins, -0.05, 0.05, nbins, -5,  5 },
    {"Zsv", /*KF*/"fZ",    true,  /*MC*/"fDecayZ",  /*KF*/"fErrZ",    "Z_{SV}", "cm",    nbins, -20,  20,  nbins, -0.05, 0.05, nbins, -5,  5 },
    {"Xpv", /*Li*/"fPosX", false, /*MC*/"fEventX",  /*KF*/"fErrPVX",  "X_{PV}", "cm",    nbins, -0.3, 0.3, nbins, -0.05, 0.05, nbins, -5,  5 },
    {"Ypv", /*Li*/"fPosY", false, /*MC*/"fEventY",  /*KF*/"fErrPVY",  "Y_{PV}", "cm",    nbins, -0.3, 0.3, nbins, -0.05, 0.05, nbins, -5,  5 },
    {"Zpv", /*Li*/"fPosZ", false, /*MC*/"fEventZ",  /*KF*/"fErrPVZ",  "Z_{PV}", "cm",    nbins, -20,  20,  nbins, -0.05, 0.05, nbins, -5,  5 },
    {"L",   /*KF*/"fL",    true,  /*MC*/"fDecayL",  /*KF*/"fDeltaL",  "L",      "cm",    nbins, -0.1, 0.2, nbins, -0.05, 0.05, nbins, -5,  5 },
    {"T",   /*KF*/"fT",    true,  /*MC*/"fDecayT",  /*KF*/"fDeltaT",  "T",      "ps",    nbins, -1,   2,   nbins, -2, 2,       nbins, -10, 10},
  };

  std::vector<std::vector<TH1D*>> hmc;  // only those mc which are matched to rec, i.e. were reconstructed
  std::vector<std::vector<TH1D*>> hrec; // only those rec which are matched to mc, i.e. correspond to signal
  std::vector<std::vector<TH1D*>> hres;
  std::vector<std::vector<TH2D*>> hcorr;
  std::vector<std::vector<TH1D*>> hpull;
  Vector2DResize(hmc, vars.size(), kNumberOfTypes);
  Vector2DResize(hrec, vars.size(), kNumberOfTypes);
  Vector2DResize(hres, vars.size(), kNumberOfTypes);
  Vector2DResize(hcorr, vars.size(), kNumberOfTypes);
  Vector2DResize(hpull, vars.size(), kNumberOfTypes);

  for(int iVar=0; iVar<vars.size(); iVar++) {
    for(int iSBType=0; iSBType<kNumberOfTypes; iSBType++) {
      hmc.at(iVar).at(iSBType) = new TH1D(("hMc_" + vars.at(iVar).name_ + "_" + dirTypes.at(iSBType)).c_str(), ("hMc_" + vars.at(iVar).name_).c_str(), vars.at(iVar).nbins_plain_, vars.at(iVar).xlow_plain_, vars.at(iVar).xup_plain_);
      hrec.at(iVar).at(iSBType) = new TH1D(("hRec_" + vars.at(iVar).name_ + "_" + dirTypes.at(iSBType)).c_str(), ("hRec_" + vars.at(iVar).name_).c_str(), vars.at(iVar).nbins_plain_, vars.at(iVar).xlow_plain_, vars.at(iVar).xup_plain_);
      hres.at(iVar).at(iSBType) = new TH1D(("hRes_" + vars.at(iVar).name_ + "_" + dirTypes.at(iSBType)).c_str(), ("hRes_" + vars.at(iVar).name_).c_str(), vars.at(iVar).nbins_res_, vars.at(iVar).xlow_res_, vars.at(iVar).xup_res_);
      hcorr.at(iVar).at(iSBType) = new TH2D(("hCorr_" + vars.at(iVar).name_ + "_" + dirTypes.at(iSBType)).c_str(), ("hCorr_" + vars.at(iVar).name_).c_str(), vars.at(iVar).nbins_plain_, vars.at(iVar).xlow_plain_, vars.at(iVar).xup_plain_, vars.at(iVar).nbins_plain_, vars.at(iVar).xlow_plain_, vars.at(iVar).xup_plain_);
      hpull.at(iVar).at(iSBType) = new TH1D(("hPull_" + vars.at(iVar).name_ + "_" + dirTypes.at(iSBType)).c_str(), ("hPull_" + vars.at(iVar).name_).c_str(), vars.at(iVar).nbins_pull_, vars.at(iVar).xlow_pull_, vars.at(iVar).xup_pull_);

      hmc.at(iVar).at(iSBType)->GetXaxis()->SetTitle((vars.at(iVar).title_ + "^{mc}, " + vars.at(iVar).unit_).c_str());
      hmc.at(iVar).at(iSBType)->GetYaxis()->SetTitle("Entries");

      hrec.at(iVar).at(iSBType)->GetXaxis()->SetTitle((vars.at(iVar).title_ + "^{rec}, " + vars.at(iVar).unit_).c_str());
      hrec.at(iVar).at(iSBType)->GetYaxis()->SetTitle("Entries");

      hres.at(iVar).at(iSBType)->GetXaxis()->SetTitle((vars.at(iVar).title_ + "^{rec} - " + vars.at(iVar).title_ + "^{mc}, " + vars.at(iVar).unit_).c_str());
      hres.at(iVar).at(iSBType)->GetYaxis()->SetTitle("Entries");

      hcorr.at(iVar).at(iSBType)->GetXaxis()->SetTitle((vars.at(iVar).title_ + "^{mc}, " + vars.at(iVar).unit_).c_str());
      hcorr.at(iVar).at(iSBType)->GetYaxis()->SetTitle((vars.at(iVar).title_ + "^{rec}, " + vars.at(iVar).unit_).c_str());
      hcorr.at(iVar).at(iSBType)->GetZaxis()->SetTitle("Entries");

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

    TTree* treeLite = fileIn->Get<TTree>((dirname + "/O2hfcandlclite").c_str());
    TTree* treeKF = fileIn->Get<TTree>((dirname + "/O2hfcandlckf").c_str());
    TTree* treeMC = fileIn->Get<TTree>((dirname + "/O2hfcandlcmc").c_str());
    treeKF->SetBranchAddress("fSigBgStatus", &sb_status);
    treeKF->SetBranchAddress("fIsSelected", &is_selected);
    for(int iVar=0; iVar<vars.size(); iVar++) {
      auto treeRec = vars.at(iVar).rec_from_kf_ ? treeKF : treeLite;
      treeRec->SetBranchAddress(vars.at(iVar).name_in_tree_rec_.c_str(), &value_rec.at(iVar));
      treeKF->SetBranchAddress(vars.at(iVar).name_in_tree_error_.c_str(), &value_err.at(iVar));
      treeMC->SetBranchAddress(vars.at(iVar).name_in_tree_mc_.c_str(), &value_mc.at(iVar));
    }

    const int Nentries = treeKF->GetEntries();

    std::cout << dirname << ", Nentries = " << Nentries << "\n";

    for(int iEntry=0; iEntry<Nentries; iEntry++) {
      treeLite->GetEntry(iEntry);
      treeKF->GetEntry(iEntry);
      treeMC->GetEntry(iEntry);

      if(is_selected<selectionFlag) continue;

      hSBType->Fill(sb_status);
      const int sb_histotype = SB_Statys2Type(sb_status);

      for(int iVar=0; iVar<vars.size(); iVar++) {
        hmc.at(iVar).at(sb_histotype)->Fill(value_mc.at(iVar));
        hrec.at(iVar).at(sb_histotype)->Fill(value_rec.at(iVar));
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
      hmc.at(iVar).at(kType)->Write();
      CD(fileOut, (histTypes.at(1) + "/" + dirTypes.at(kType)).c_str());
      hrec.at(iVar).at(kType)->Write();
      CD(fileOut, (histTypes.at(2) + "/" + dirTypes.at(kType)).c_str());
      hres.at(iVar).at(kType)->Write();
      CD(fileOut, (histTypes.at(3) + "/" + dirTypes.at(kType)).c_str());
      hcorr.at(iVar).at(kType)->Write();
      CD(fileOut, (histTypes.at(4) + "/" + dirTypes.at(kType)).c_str());
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
