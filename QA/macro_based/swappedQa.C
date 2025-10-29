std::string ReadNthLine(const std::string& fileName);
std::vector<std::string> GetDFNames(const std::string& fileName);

/// Channels taken from here: https://github.com/AliceO2Group/O2Physics/blob/87be5da87be8bcef56dccf64b5d960e7f1b7545d/PWGHF/Core/DecayChannels.h#L61-L95
/// @brief 3-prong candidates: main channels
enum DecayChannelMain : int8_t {
  // D+
  DplusToPiKPi = 1,    // π+ K− π+
  DplusToPiKPiPi0 = 2, // π+ K− π+ π0
  DplusToPiPiPi = 3,   // π+ π− π+
  DplusToPiKK = 4,     // π+ K− K+
  // Ds+
  DsToPiKK = 5,      // π+ K− K+
  DsToPiKKPi0 = 6,   // π+ K− K+ π0
  DsToPiPiK = 7,     // π+ π− K+
  DsToPiPiPi = 8,    // π+ π− π+
  DsToPiPiPiPi0 = 9, // π+ π− π+ π0
  // D*+
/*[x]*/  DstarToPiKPi = 10,       // π+ K− π+ (from [(D0 → π+ K−) π+])
/*[x]*/  DstarToPiKPiPi0 = 11,    // π+ K− π+ π0
/*[ ]*/  DstarToPiKPiPi0Pi0 = 12, // π+ K− π+ π0 π0
/*[x]*/  DstarToPiKK = 13,        // π+ K− K+
/*[ ]*/  DstarToPiKKPi0 = 14,     // π+ K− K+ π0
/*[x]*/  DstarToPiPiPi = 15,      // π+ π− π+
/*[x]*/  DstarToPiPiPiPi0 = 16,   // π+ π− π+ π0
  // Λc+
  LcToPKPi = 17,    // p K− π+
  LcToPKPiPi0 = 18, // p K− π+ π0
  LcToPPiPi = 19,   // p π− π+
  LcToPKK = 20,     // p K− K+
  // Ξc+
  XicToPKPi = 21,  // p K− π+
  XicToPKK = 22,   // p K− K+
  XicToSPiPi = 23, // Σ+ π− π+
  //
  NChannelsMain = XicToSPiPi // last channel
};

struct Mother {
  int id_;
  std::string name_;
  std::string greek_name_;
};

enum MotherParticle : int {
  Dplus = 0,
  Ds,
  Dstar,
  Lc,
  Xic,
  NMotherParticles
};

std::vector<Mother> Mothers {
  {Dplus, "Dplus", "D^{+}"          },
  {Ds,    "Ds",    "D_{s}^{+}"      },
  {Dstar, "Dstar", "D^{*+}"         },
  {Lc,    "Lc",    "#Lambda_{c}^{+}"},
  {Xic,   "Xic",   "#Xi_{c}^{+}"    }
};

struct Decay {
  int id_;
  Mother mother_;
  std::string daughters_;
  std::string greek_daughters_;
};
std::vector<Decay> Decays {
  {DplusToPiKPi,       Mothers[Dplus], "PiKPi",       "#pi^{+}K^{#minus}#pi^{+}"              },
  {DplusToPiKPiPi0,    Mothers[Dplus], "PiKPiPi0",    "#pi^{+}K^{#minus}#pi^{+}#pi^{0}"       },
  {DplusToPiPiPi,      Mothers[Dplus], "PiPiPi",      "#pi^{+}#pi^{#minus}#pi^{+}"            },
  {DplusToPiKK,        Mothers[Dplus], "PiKK",        "#pi^{+}K^{#minus}K^{+}"                },

  {DsToPiKK,           Mothers[Ds],    "PiKK",        "#pi^{+}K^{#minus}K^{+}"                },
  {DsToPiKKPi0,        Mothers[Ds],    "PiKKPi0",     "#pi^{+}K^{#minus}K^{+}#pi^{0}"         },
  {DsToPiPiK,          Mothers[Ds],    "PiPiK",       "#pi^{+}#pi^{#minus}K^{+}"              },
  {DsToPiPiPi,         Mothers[Ds],    "PiPiPi",      "#pi^{+}#pi^{#minus}#pi^{+}"            },
  {DsToPiPiPiPi0,      Mothers[Ds],    "PiPiPiPi0",   "#pi^{+}#pi^{#minus}#pi^{+}#pi^{0}"     },

  {DstarToPiKPi,       Mothers[Dstar], "PiKPi",       "#pi^{+}K^{#minus}#pi^{+}"              },
  {DstarToPiKPiPi0,    Mothers[Dstar], "PiKPiPi0",    "#pi^{+}K^{#minus}#pi^{+}#pi^{0}"       },
  {DstarToPiKPiPi0Pi0, Mothers[Dstar], "PiKPiPi0Pi0", "#pi^{+}K^{#minus}#pi^{+}#pi^{0}#pi^{0}"},
  {DstarToPiKK,        Mothers[Dstar], "PiKK",        "#pi^{+}K^{#minus}K^{+}"                },
  {DstarToPiKKPi0,     Mothers[Dstar], "PiKKPi0",     "#pi^{+}K^{#minus}K^{+}#pi^{0}"         },
  {DstarToPiPiPi,      Mothers[Dstar], "PiPiPi",      "#pi^{+}#pi^{#minus}#pi^{+}"            },
  {DstarToPiPiPiPi0,   Mothers[Dstar], "PiPiPiPi0",   "#pi^{+}#pi^{#minus}#pi^{+}#pi^{0}"     },

  {LcToPKPi,           Mothers[Lc],    "PKPi",        "pK^{#minus}#pi^{+}"                    },
  {LcToPKPiPi0,        Mothers[Lc],    "PKPiPi0",     "pK^{#minus}#pi^{+}#pi^{0}"             },
  {LcToPPiPi,          Mothers[Lc],    "PPiPi",       "p#pi^{#minus}#pi^{+}"                  },
  {LcToPKK,            Mothers[Lc],    "PKK",         "pK^{#minus}K^{+}"                      },

  {XicToPKPi,          Mothers[Xic],   "PKPi",        "pK^{#minus}#pi^{+}"                    },
  {XicToPKK,           Mothers[Xic],   "PKK",         "pK^{#minus}K^{+}"                      },
  {XicToSPiPi,         Mothers[Xic],   "SPiPi",       "#Sigma^{+}#pi^{#minus}#pi^{+}"         },
};

void swappedQa(const std::string& filenameIn, const int nFiles) {
  const size_t nDecays{Decays.size()};
  std::vector<TH1*> histos;
  histos.resize(nDecays);
  std::vector<std::string> cuts;
  cuts.resize(nDecays);
  std::vector<std::string> histoNames;
  histoNames.resize(nDecays);

  for(size_t iDecay=0; iDecay<nDecays; ++iDecay) {
    const std::string histoName = Decays.at(iDecay).mother_.name_ + "To" + Decays.at(iDecay).daughters_;
    histoNames.at(iDecay) = histoName.c_str();
    histos.at(iDecay) = new TH1D(histoName.c_str(), (Decays.at(iDecay).mother_.greek_name_ + "#rightarrow" + Decays.at(iDecay).greek_daughters_).c_str(), 4, -2, 2);
    histos.at(iDecay)->GetXaxis()->SetTitle("isCandidateSwapped");
    histos.at(iDecay)->GetYaxis()->SetTitle("Entries");
    cuts.at(iDecay) = "fFlagMc == " + std::to_string(-Decays.at(iDecay).id_) + " || fFlagMc == " + std::to_string(Decays.at(iDecay).id_);
  }

  for(int iFile=1; iFile<=nFiles; ++ iFile) {
    std::cout << "Processing iFile = " << iFile << "\n";
    const std::string filename = ReadNthLine(filenameIn + ":" + std::to_string(iFile));
    TFile* fileIn = TFile::Open(filename.c_str(), "read");
    const std::vector<std::string> dirNames = GetDFNames(filename);

    for(const auto& dirName : dirNames) {
      std::cout << "Processing " << dirName << "\n";
      TTree* tree = fileIn->Get<TTree>((dirName + "/O2hfcandlclite").c_str());

      for(size_t iDecay=0; iDecay<nDecays; ++iDecay) {
        tree->Draw(Form("fIsCandidateSwapped>>+%s", histoNames.at(iDecay).c_str()), cuts.at(iDecay).c_str(), "goff");
      } // nDecays
    } // dirNames
    std::cout << "\n";
  } // nFiles


  TFile* fileOut = TFile::Open("swappedQa.root", "recreate");
  for(size_t iDecay=0; iDecay<nDecays; ++iDecay) {
    histos.at(iDecay)->Write();
  }
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
    if(dirname.substr(0, 2) != "DF") continue;
    result.emplace_back(dirname);
  }
  fileIn->Close();

  return result;
}

std::string ReadNthLine(const std::string& fileName) {
  if(fileName.find(':') == std::string::npos) return fileName;

  std::string result;
  const size_t colonPosition = fileName.find(':');
  const std::string fileListName = fileName.substr(0, colonPosition);
  const std::string fileLineNumberStr = fileName.substr(colonPosition + 1);
  const int fileLineNumberInt = std::stoi(fileLineNumberStr);

  std::ifstream fileList(fileListName);
  if (!fileList.is_open()) throw std::runtime_error("ReadNthLine() - the fileList " + fileListName + " is missing!");

  for(size_t iLine=0; iLine<fileLineNumberInt; ++iLine) {
    if (!std::getline(fileList, result)) throw std::runtime_error("ReadNthLine() - the EOF of fileList " + fileListName + " reached before line " + fileLineNumberStr);
  }

  return result;
}
