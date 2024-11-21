int SB_Statys2Type(int status);
void CD(TFile* file, std::string dir);
enum eType : short {
  kPrompt = 0,
  kNonPrompt,
  kBackground,
  kWrongSwap,
  kImpossible,
  kNumberOfTypes
};

void treeKF_qa(const std::string& fileName, int selectionFlag=1) {

  TFile* fileIn = TFile::Open(fileName.c_str());
  if(fileIn == nullptr) {
    throw std::runtime_error("fileIn == nullptr");
  }

  TFile* fileOut = TFile::Open("treeKF_qa.root", "recreate");

  std::vector<std::string> dirTypes{"prompt", "nonprompt", "background", "wrongswap", "impossible"};

  TH1D* hSBType = new TH1D("hSBType", "hSBType", 10, -2.5, 7.5);

  struct Variable {
    std::string name_;
    std::string name_in_tree_;
    std::string xaxis_title_;
    int nbins_;
    float xlow_;
    float xup_;
  };

  std::vector<Variable> vars {
    {"Mass",         "fMassInv",           "m_{pK#pi}, GeV/c^{2}", 600, 1.98, 2.58},
    {"Pt",           "fPt",                "p_{T}, GeV/c",         600, 0,    12  },
    {"Chi2prim_p",   "fChi2PrimProton",    "#chi^{2}_{prim}{p}",   100, 0,    60  },
    {"Chi2prim_K",   "fChi2PrimKaon",      "#chi^{2}_{prim}{K}",   100, 0,    60  },
    {"Chi2prim_pi",  "fChi2PrimPion",      "#chi^{2}_{prim}{#pi}", 100, 0,    60  },
    {"Chi2geo_p_pi", "fChi2geoProtonPion", "#chi^{2}_{geo}{p#pi}", 100, 0,    4   },
    {"Chi2geo_p_K",  "fChi2geoProtonKaon", "#chi^{2}_{geo}{pK}",   100, 0,    4   },
    {"Chi2geo_K_pi", "fChi2geoPionKaon",   "#chi^{2}_{geo}{K#pi}", 100, 0,    4   },
    {"DCA_p_pi",     "fDCAProtonPion",     "DCA{p#pi}, cm",        100, 0,    0.1 },
    {"DCA_p_K",      "fDCAProtonKaon",     "DCA{pK}, cm",          100, 0,    0.1 },
    {"DCA_K_pi",     "fDCAPionKaon",       "DCA{K#pi}, cm",        100, 0,    0.1 },
    {"Chi2geo",      "fChi2geo",           "#chi^{2}_{geo}",       100, 0,    10  },
    {"Chi2topo",     "fChi2topo",          "#chi^{2}_{topo}",      100, 0,    20  },
    {"LdL",          "fLdL",               "L/#Delta L",           100, 0,    10  },
    {"L",            "fL",                 "L, cm",                100, 0,    0.5 },
    {"T",            "fT",                 "T, ps",                400, 0,    5   },
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

  std::vector<float> value;
  value.resize(vars.size());
  int sb_status, is_selected;

  auto lok = fileIn->GetListOfKeys();
  fileIn->cd();
  for(const auto& k : *lok) {
    const std::string dirname = k->GetName();
    if(dirname == "parentFiles") continue;

    TTree* treeKF = fileIn->Get<TTree>((dirname + "/O2hfcandlckf").c_str());
    treeKF->SetBranchAddress("fSigBgStatus", &sb_status);
    treeKF->SetBranchAddress("fIsSelected", &is_selected);
    for(int iVar=0; iVar<vars.size(); iVar++) {
      treeKF->SetBranchAddress(vars.at(iVar).name_in_tree_.c_str(), &value.at(iVar));
    }

    const int Nentries = treeKF->GetEntries();
    for(int iEntry=0; iEntry<Nentries; iEntry++) {
      treeKF->GetEntry(iEntry);

      if(is_selected<selectionFlag) continue;

      hSBType->Fill(sb_status);
      const int sb_histotype = SB_Statys2Type(sb_status);

      for(int iVar=0; iVar<vars.size(); iVar++) {
        hvar.at(iVar).at(sb_histotype)->Fill(value.at(iVar));
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
  else return kImpossible;
}

void CD(TFile* file, std::string dir) {
  if(file == nullptr) throw std::runtime_error("CD() - file is nullptr");

  if(file->GetDirectory(dir.c_str()) == nullptr) file->mkdir(dir.c_str());
  file->cd(dir.c_str());
}
