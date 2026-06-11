void lumiRead(const std::string& fileList) {
  ifstream in(fileList.c_str());
  if (!in.is_open()) {
    std::cout << "Cannot open file list!" << endl;
    return;
  }

  std::string filename;

  double counterTVX{};
  double counterTCE{};
  double counterTVXafterBCcuts{};
  double counterTCEafterBCcuts{};
  double lumiTVX{};
  double lumiTCE{};
  double lumiTVXafterBCcuts{};
  double lumiTCEafterBCcuts{};

  int nFile{};

  while (getline(in, filename)) {
    if(nFile%10 == 0) std::cout << "File no. " << nFile << "\n";

    if (filename.empty()) continue;

    TFile *f = TFile::Open(filename.c_str(), "READ");

    if (!f || f->IsZombie()) {
      std::cout << "Error opening file: " << filename << endl;
      continue;
    }

    TH1* hCounterTVX = f->Get<TH1>("eventselection-run3/luminosity/hCounterTVX");
    TH1* hCounterTCE = f->Get<TH1>("eventselection-run3/luminosity/hCounterTCE");
    TH1* hCounterTVXafterBCcuts = f->Get<TH1>("eventselection-run3/luminosity/hCounterTVXafterBCcuts");
    TH1* hCounterTCEafterBCcuts = f->Get<TH1>("eventselection-run3/luminosity/hCounterTCEafterBCcuts");
    if(hCounterTVX == nullptr || hCounterTCE == nullptr || hCounterTVXafterBCcuts == nullptr || hCounterTCEafterBCcuts == nullptr) throw std::runtime_error("At least one of the histograms is nullptr");
    TH1* hLumiTVX = f->Get<TH1>("eventselection-run3/luminosity/hLumiTVX");
    TH1* hLumiTCE = f->Get<TH1>("eventselection-run3/luminosity/hLumiTCE");
    TH1* hLumiTVXafterBCcuts = f->Get<TH1>("eventselection-run3/luminosity/hLumiTVXafterBCcuts");
    TH1* hLumiTCEafterBCcuts = f->Get<TH1>("eventselection-run3/luminosity/hLumiTCEafterBCcuts");
    if(hLumiTVX == nullptr || hLumiTCE == nullptr || hLumiTVXafterBCcuts == nullptr || hLumiTCEafterBCcuts == nullptr) throw std::runtime_error("At least one of the histograms is nullptr");

    counterTVX += hCounterTVX->GetEntries();
    counterTCE += hCounterTCE->GetEntries();
    counterTVXafterBCcuts += hCounterTVXafterBCcuts->GetEntries();
    counterTCEafterBCcuts += hCounterTCEafterBCcuts->GetEntries();

    lumiTVX += hLumiTVX->Integral();
    lumiTCE += hLumiTCE->Integral();
    lumiTVXafterBCcuts += hLumiTVXafterBCcuts->Integral();
    lumiTCEafterBCcuts += hLumiTCEafterBCcuts->Integral();

    f->Close();
    delete f;

    ++nFile;
  }

  std::cout << "\n";
  std::cout << "counterTVX = " << counterTVX << "\n";
  std::cout << "counterTCE = " << counterTCE << "\n";
  std::cout << "counterTVXafterBCcuts = " << counterTVXafterBCcuts << "\n";
  std::cout << "counterTCEafterBCcuts = " << counterTCEafterBCcuts << "\n";
  std::cout << "lumiTVX = " << lumiTVX << "\n";
  std::cout << "lumiTCE = " << lumiTCE << "\n";
  std::cout << "lumiTVXafterBCcuts = " << lumiTVXafterBCcuts << "\n";
  std::cout << "lumiTCEafterBCcuts = " << lumiTCEafterBCcuts << "\n";

  in.close();
}
