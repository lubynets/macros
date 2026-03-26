TH1* EvaluateEfficiency(const TH1* histoNum, const TH1* histoDen) {
  TH1* histoEff = dynamic_cast<TH1*>(histoNum->Clone());
  histoEff->Sumw2();
  histoEff->Divide(histoDen);

  return histoEff;
}

void fitEff() {
  TFile* fileYield = TFile::Open("ct_yield_qa.HF_LHC24h1b_All.595984.root", "read");
  if(fileYield == nullptr) throw std::runtime_error("fileYield == nullptr");
  TH1* hYieldGen = fileYield->Get<TH1>("prompt/hGen");
  TH1* hYieldCand = fileYield->Get<TH1>("prompt/hCand");
  TH1* hYieldSim = fileYield->Get<TH1>("prompt/hSim");
  if(hYieldGen == nullptr || hYieldCand == nullptr || hYieldSim == nullptr) throw std::runtime_error("hYieldGen == nullptr || hYieldCand == nullptr || hYieldSim == nullptr");

  TH1* hEffSim = EvaluateEfficiency(hYieldSim, hYieldGen);
  TH1* hEffCand = EvaluateEfficiency(hYieldCand, hYieldGen);

  auto stitched = [](double *x, double *p) {
      double xx = x[0];

      double a = p[0];
      double x0 = p[1];
      double sigma = p[2];
      double A = p[3];
      double m = p[4];

      double erf_at_a = 0.5 * A * (1 + TMath::Erf((a - x0)/sigma));
      double b = erf_at_a - m * a;

      if (xx < a) {
          return 0.5 * A * (1 + TMath::Erf((xx - x0)/sigma));
      } else {
          return m * xx + b;
      }
  };

  TFile* fileOut = TFile::Open("eff_fit.root", "recreate");

  TF1 *f = new TF1("f", stitched, 0., 2., 5);
  f->FixParameter(0, 0.6);

  hEffSim->Fit(f, "R", "");
  f->Write("effSim");

  hEffCand->Fit(f, "R", "");
  f->Write("effCand");

  fileOut->Close();
  fileYield->Close();
}
