void pid_qa(const char* inFileName) {

  std::ifstream inputList(inFileName);

  if (!inputList.is_open()) {
    Error("pid_qa", "Cannot open input list %s", inFileName);
    return;
  }

  TChain chain;

  std::string fileName;

  while (std::getline(inputList, fileName)) {

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

      TString treePath = dirName + "/O2qapidtpc";

      chain.Add(Form("%s/%s", fileName.c_str(), treePath.Data()));
    }
  }

  ROOT::RDataFrame df(chain);

  struct Species {
    int pid;
    const char* name;
    const char* symbol;
  };

  const std::vector<Species> species = {
    {0, "El", "e"},
    {2, "Pi", "#pi"},
    {3, "Ka", "K"},
    {4, "Pr", "p"}
  };

  struct HistDef {
    const char* branch;
    const char* directory;
    int nBinsY;
    double yMin;
    double yMax;
    const char* title;
  };

  const std::vector<HistDef> histDefs = {
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

  const int nBinsP = 3000;
  const double pMin = 0.01;
  const double pMax = 20.;

  std::vector<double> pBins(nBinsP + 1);

  const double logMin = std::log10(pMin);
  const double logMax = std::log10(pMax);

  for (int i = 0; i <= nBinsP; ++i) {
    pBins[i] = std::pow(10., logMin + (logMax - logMin) * i / nBinsP);
  }

  std::map<std::string, ROOT::RDF::RResultPtr<TH1>> histos;

  for (const auto& hd : histDefs) {

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

  for (const auto& hd : histDefs) {

    output.mkdir(hd.directory);
    output.cd(hd.directory);

    for (const auto& s : species) {
      histos[static_cast<std::string>(hd.directory) + "_" + static_cast<std::string>(s.name)]->Write(s.name);
    }

    output.cd();
  }

  output.Close();
}
