void pid_qa(const char* inFileName) {

  std::string inputName(inFileName);
  int selectedLine = 0;  // 0 means process all lines

  auto colonPos = inputName.rfind(':');

  if (colonPos != std::string::npos) {
    std::string lineStr = inputName.substr(colonPos + 1);

    if (!lineStr.empty() &&
        std::all_of(lineStr.begin(), lineStr.end(), ::isdigit)) {
      selectedLine = std::stoi(lineStr);
      inputName = inputName.substr(0, colonPos);
    }
  }

  std::ifstream inputList(inputName);

  if (!inputList.is_open()) {
    Error("pid_qa", "Cannot open input list %s", inputName.c_str());
    return;
  }

  TChain chain;

  std::string fileName;

  int iFile{0};
  while (std::getline(inputList, fileName)) {
  ++iFile;

  if (selectedLine > 0 && iFile != selectedLine)
    continue;

  std::cout << "Capturing file no. " << iFile << std::endl;

  if (fileName.empty())
    continue;

    TFile input(fileName.c_str());

    if (input.IsZombie())
      continue;

    TIter next(input.GetListOfKeys());
    TKey* key;

    while ((key = (TKey*)next())) {

      if (strcmp(key->GetClassName(), "TDirectoryFile") != 0)
        continue;

      TString dirName = key->GetName();

      if (!dirName.BeginsWith("DF_"))
        continue;

      TString treePath = dirName + "/O2pidtpcdiag";

      chain.Add(Form("%s/%s", fileName.c_str(), treePath.Data()));
    }
  }

  // ============================================================
  const double cutY = 70.;
  const double cutX = 0.4;
  const double leftEdgeX = 0.1;

  ROOT::RDataFrame df0(chain);

  auto df = df0.Filter(
    [cutX, cutY, leftEdgeX](float p, float sigma) {
      return p > 0. &&
            sigma < cutY *
            std::log10(cutX / p) / std::log10(cutX / leftEdgeX);
    },
    {"fTPCInnerParam", "fExpSigma"}
  );
  // ============================================================

  struct Species {
    int pid;
    const char* name;
    const char* symbol;
  };

  const std::vector<Species> species = {
//     {0, "El", "e"},
//     {2, "Pi", "#pi"},
//     {3, "Ka", "K"},
    {4, "Pr", "p"}
  };

  struct Hist2DDef {
    const char* branch;
    const char* directory;
    int nBinsY;
    double yMin;
    double yMax;
    const char* title;
  };

  const std::vector<Hist2DDef> hist2DDefs = {
    {
      "fNSigmaTpc",
      "nsigma",
      401, -10.025, 10.025,
      "N_{#sigma}^{TPC}(%s)"
    },
    {
      "fDedxExpected",
      "expected",
      1000, 0., 1000.,
      "d#it{E}/d#it{x}(%s) A.U."
    },
    {
      "fDedxDiff",
      "delta",
      200, -1000., 1000.,
      "d#it{E}/d#it{x} - d#it{E}/d#it{x}(%s)"
    },
    {
      "fExpSigma",
      "expsigma",
      200, 0., 200.,
      "Exp_{#sigma}^{TPC}(%s)"
    }
  };

  struct Hist1DDef {
    const char* branch;
    int nBins;
    double xMin;
    double xMax;
  };

  const std::vector<Hist1DDef> hist1DDefs = {
    {"fFt0Occ", 200, 0., 100000.},
    {"fHadronicRate", 200, 15., 25.},
    {"fNormMultTPC", 200, 0., 2.},
    {"fNormNClustersTPC", 200, 0., 2.},
    {"fPhi", 200, 0., 2*3.1415926},
    {"fTgl", 200, -1.5, 1.5},
    {"fSigned1Pt", 200, -10., 10.},
    {"fY", 200, -1.2, 1.2},
    {"fP", 200, 0., 100.},
    {"fTPCInnerParam", 200, 0., 100.}
  };

  const int nBinsP = 2000;
  const double pMin = 0.1;
  const double pMax = 20.;

  std::vector<double> pBins(nBinsP + 1);

  const double logMin = std::log10(pMin);
  const double logMax = std::log10(pMax);

  for (int i = 0; i <= nBinsP; ++i) {
    pBins[i] = std::pow(10., logMin + (logMax - logMin) * i / nBinsP);
  }

  std::map<std::string, ROOT::RDF::RResultPtr<TH1>> histos;

  std::map<std::string, ROOT::RDF::RResultPtr<TH1D>> histos1D;

  for (const auto& hd : hist1DDefs) {

    auto h = df.Histo1D(
      {
        hd.branch,
        hd.branch,
        hd.nBins,
        hd.xMin,
        hd.xMax
      },
      hd.branch
    );

    h->GetXaxis()->SetTitle(hd.branch);

    histos1D[hd.branch] = h;
  }

  for (const auto& hd : hist2DDefs) {

    for (const auto& s : species) {

      TString histName = Form("%s_%s", hd.directory, s.name);

      auto h =
        df.Filter(Form("fPidIndex == %d", s.pid))
          .Histo2D(
            {
              histName.Data(),
              Form(hd.title, s.symbol),
              nBinsP, pBins.data(),
              hd.nBinsY, hd.yMin, hd.yMax
            },
            std::strcmp(hd.branch, "fNSigmaTpc") == 0 ? "fP" : "fTPCInnerParam",
            hd.branch);

      h->GetXaxis()->SetTitle("#it{p}/|Z| (GeV/#it{c})");
      h->GetYaxis()->SetTitle(h->GetTitle());

      histos[histName.Data()] = h;

    }
  }

  TFile output("pid_qa.root", "recreate");

  for (const auto& hd : hist2DDefs) {

    output.mkdir(hd.directory);
    output.cd(hd.directory);

    for (const auto& s : species) {
      histos[static_cast<std::string>(hd.directory) + "_" + static_cast<std::string>(s.name)]->Write(s.name);
    }

    output.cd();
  }

  output.mkdir("trainingVars");
  output.cd("trainingVars");

  for (const auto& hd : hist1DDefs) {
    histos1D[hd.branch]->Write();
  }

  output.cd();

  output.Close();
}
