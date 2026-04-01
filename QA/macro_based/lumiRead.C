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

    counterTVX += hCounterTVX->GetEntries();
    counterTCE += hCounterTCE->GetEntries();
    counterTVXafterBCcuts += hCounterTVXafterBCcuts->GetEntries();
    counterTCEafterBCcuts += hCounterTCEafterBCcuts->GetEntries();

    f->Close();
    delete f;

    ++nFile;
  }

  std::cout << "\n";
  std::cout << "counterTVX = " << counterTVX << "\n";
  std::cout << "counterTCE = " << counterTCE << "\n";
  std::cout << "counterTVXafterBCcuts = " << counterTVXafterBCcuts << "\n";
  std::cout << "counterTCEafterBCcuts = " << counterTCEafterBCcuts << "\n";

  in.close();
}
