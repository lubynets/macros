int SB_Statys2Type(int status);
void CD(TFile* file, std::string dir);
enum eType : short {
  kPrompt = 0,
  kNonPrompt,
  kBackground,
  kWrongSwap,
  kImpossible,
  kData,
  kNumberOfTypes
};

void treeKF_qa(const std::string& fileName, int selectionFlag=1) {

  TFile* fileIn = TFile::Open(fileName.c_str());
  if(fileIn == nullptr) {
    throw std::runtime_error("fileIn == nullptr");
  }

  TFile* fileOut = TFile::Open("treeKF_qa.root", "recreate");

  std::vector<std::string> dirTypes{"prompt", "nonprompt", "background", "wrongswap", "impossible", "data"};

  TH1D* hSBType = new TH1D("hSBType", "hSBType", 8, -2.5, 5.5);
  const char* sb_titles[8] = {"", "Impossible", "Background", "Prompt", "Non-prompt", "Wrong swap", "Data", ""};
  for (int i=1;i<=8;i++) {
    hSBType->GetXaxis()->SetBinLabel(i,sb_titles[i-1]);
  }

  struct Variable {
    std::string name_;
    std::string name_in_tree_;
    bool is_from_kf_;
    bool is_float_;
    std::string xaxis_title_;
    int nbins_;
    float xlow_;
    float xup_;
  };

  std::vector<Variable> vars {
    {"Mass",         "fMassInv",               true,  true,  "m_{pK#pi}, GeV/c^{2}", 600, 1.98, 2.58},
    {"Pt",           "fPt",                    true,  true,  "p_{T}, GeV/c",         600, 0,    12  },
    {"Chi2prim_p",   "fChi2PrimProton",        true,  true,  "#chi^{2}_{prim}{p}",   100, 0,    60  },
    {"Chi2prim_K",   "fChi2PrimKaon",          true,  true,  "#chi^{2}_{prim}{K}",   100, 0,    60  },
    {"Chi2prim_pi",  "fChi2PrimPion",          true,  true,  "#chi^{2}_{prim}{#pi}", 100, 0,    60  },
    {"Chi2geo_p_pi", "fChi2geoProtonPion",     true,  true,  "#chi^{2}_{geo}{p#pi}", 100, 0,    4   },
    {"Chi2geo_p_K",  "fChi2geoProtonKaon",     true,  true,  "#chi^{2}_{geo}{pK}",   100, 0,    4   },
    {"Chi2geo_K_pi", "fChi2geoPionKaon",       true,  true,  "#chi^{2}_{geo}{K#pi}", 100, 0,    4   },
    {"DCA_p_pi",     "fDCAProtonPion",         true,  true,  "DCA{p#pi}, cm",        100, 0,    0.1 },
    {"DCA_p_K",      "fDCAProtonKaon",         true,  true,  "DCA{pK}, cm",          100, 0,    0.1 },
    {"DCA_K_pi",     "fDCAPionKaon",           true,  true,  "DCA{K#pi}, cm",        100, 0,    0.1 },
    {"Chi2geo",      "fChi2geo",               true,  true,  "#chi^{2}_{geo}",       100, 0,    10  },
    {"Chi2topo",     "fChi2topo",              true,  true,  "#chi^{2}_{topo}",      100, 0,    20  },
    {"LdL",          "fLdL",                   true,  true,  "L/#Delta L",           100, 0,    10  },
    {"L",            "fL",                     true,  true,  "L, cm",                100, 0,    0.5 },
    {"T",            "fT",                     true,  true,  "T, ps",                400, 0,    5   },
    {"nPCPV",        "fNProngsContributorsPV", false, false, "nPCPV",                6,   -1  , 5   },
  };

  std::vector<std::vector<TH1D*>> hvar;
  hvar.resize(vars.size());
  for(auto& h : hvar) {
    h.resize(kNumberOfTypes);
  }

  for(int iVar=0; iVar<vars.size(); iVar++) {
    for(int iSBType=0; iSBType<kNumberOfTypes; iSBType++) {
      hvar.at(iVar).at(iSBType) = new TH1D(("h" + vars.at(iVar).name_ + "_" + dirTypes.at(iSBType)).c_str(), ("h" + vars.at(iVar).name_).c_str(), vars.at(iVar).nbins_, vars.at(iVar).xlow_, vars.at(iVar).xup_);
      hvar.at(iVar).at(iSBType)->GetXaxis()->SetTitle(vars.at(iVar).xaxis_title_.c_str());
      hvar.at(iVar).at(iSBType)->GetYaxis()->SetTitle("Entries");
    }
  }

  std::vector<float> value_float;
  std::vector<uint8_t> value_int;
  value_float.resize(vars.size());
  value_int.resize(vars.size());
  int sb_status, is_selected;

  auto lok = fileIn->GetListOfKeys();
  fileIn->cd();
  for(const auto& k : *lok) {
    const std::string dirname = k->GetName();
    if(dirname == "parentFiles") continue;

    TTree* treeKF = fileIn->Get<TTree>((dirname + "/O2hfcandlckf").c_str());
    TTree* treeLite = fileIn->Get<TTree>((dirname + "/O2hfcandlclite").c_str());
    treeKF->SetBranchAddress("fSigBgStatus", &sb_status);
    treeKF->SetBranchAddress("fIsSelected", &is_selected);
    for(int iVar=0; iVar<vars.size(); iVar++) {
      auto treeRec = vars.at(iVar).is_from_kf_ ? treeKF : treeLite;
      if(vars.at(iVar).is_float_) treeRec->SetBranchAddress(vars.at(iVar).name_in_tree_.c_str(), &value_float.at(iVar));
      else                        treeRec->SetBranchAddress(vars.at(iVar).name_in_tree_.c_str(), &value_int.at(iVar));
    }

    const int Nentries = treeKF->GetEntries();
    for(int iEntry=0; iEntry<Nentries; iEntry++) {
      treeKF->GetEntry(iEntry);
      treeLite->GetEntry(iEntry);

      if(is_selected<selectionFlag) continue;

      hSBType->Fill(sb_status == -999 ? 4 : sb_status);
      const int sb_histotype = SB_Statys2Type(sb_status);

      for(int iVar=0; iVar<vars.size(); iVar++) {
        auto value = vars.at(iVar).is_float_ ? value_float.at(iVar) : static_cast<int>(value_int.at(iVar));
        hvar.at(iVar).at(sb_histotype)->Fill(value);
      }
    }
  }

  fileOut->cd();
  hSBType->Write();
  for(int kType=0; kType<kNumberOfTypes; kType++) {
    CD(fileOut, dirTypes.at(kType));
    for(int iVar=0; iVar<vars.size(); iVar++) {
      hvar.at(iVar).at(kType)->Write();
    }
  }

  fileOut->Close();
  fileIn->Close();
}

int SB_Statys2Type(int status) {
  if      (status == 0) return kBackground;
  else if (status == 1) return kPrompt;
  else if (status == 2) return kNonPrompt;
  else if (status == 3) return kWrongSwap;
  else if (status == -999) return kData;
  else return kImpossible;
}

void CD(TFile* file, std::string dir) {
  if(file == nullptr) throw std::runtime_error("CD() - file is nullptr");

  if(file->GetDirectory(dir.c_str()) == nullptr) file->mkdir(dir.c_str());
  file->cd(dir.c_str());
}
