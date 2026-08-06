void analyticExpoFit() {
  const std::string fileName{"/home/oleksii/alidir/working/systematics/bdtEff/standard/HF_LHC24h1b_All/667716/cutVar/keppAsIs/CutVarLc.merged.root"};
  const std::string histoName{"hCorrYieldsPrompt"};

  TFile* fileIn = TFile::Open(fileName.c_str());
  TH1* histoIn = fileIn->Get<TH1>(histoName.c_str());
  histoIn->Scale(1., "width");

  const int nBins = histoIn->GetNbinsX();

  std::vector<double> x, y, s, w, eta;
  double S{}, Sx{}, Sy{}, Sxx{}, Sxy{};

  for(int iBin=1; iBin<=nBins; ++iBin) {
    x.push_back(histoIn->GetBinCenter(iBin));
    y.push_back(std::log(histoIn->GetBinContent(iBin)));
    s.push_back(histoIn->GetBinError(iBin) / histoIn->GetBinContent(iBin));
    w.push_back(1. / s.back() / s.back());

    S += w.back();
    Sx += w.back() * x.back();
    Sy += w.back() * y.back();
    Sxx += w.back() * x.back() * x.back();
    Sxy += w.back() * x.back() * y.back();

    histoIn->SetBinContent(iBin, y.back());
    histoIn->SetBinError(iBin, s.back());

  }

  const double Delta = S*Sxx - Sx*Sx;
  const double xw = Sx / S;

  const double A = (S*Sxy - Sx*Sy) / Delta;
  const double VarA = S / Delta;
  const double sA = std::sqrt(VarA);

  std::cout << "tau = " << -1./A << " +- " << sA/A/A << "\n";

  std::vector<double> V{}, f{};
  for(int i=0; i<nBins; ++i) {
    V.push_back((S*S*(x[i] - xw)*(x[i] - xw)) / Delta / Delta / s[i] / s[i]);
    f.push_back(V.back() / VarA);

    // ----- delete point i -----

    const double Sm   = S   - w[i];
    const double Sxm  = Sx  - w[i]*x[i];
    const double Sxxm = Sxx - w[i]*x[i]*x[i];

    const double Deltam = Sm*Sxxm - Sxm*Sxm;

    const double VarAm = Sm / Deltam;
    const double sAm   = std::sqrt(VarAm);

    // relative loss of precision
    eta.push_back((sAm - sA)/sA);


    std::cout
    << i+1
    << "  x=" << x[i]
    << "  dx=" << x[i]-xw
    << "  sigma=" << s[i]
    << "  weight=" << w[i]
    << "  score="
    << (x[i]-xw)*(x[i]-xw)/s[i]/s[i]
    << "  V=" << V[i]
    << "  f=" << f[i]
    << "  eta=" << eta[i]
    << " (" << 100.*eta[i] << "%)"
    << std::endl;

  }

  histoIn->Draw();
}
