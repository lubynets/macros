const double a = 6.81, b = 59.24;
const double c = 0.082, d = 0.151;
const double sqrtSNN = 5360;
const double mass = 0.140;

double rho(double pT);

void tsallis() {
  TGraph gr;
  gr.GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
//   gr.GetYaxis()->SetTitle("tsallisCharged");
  gr.GetYaxis()->SetTitle("tsallisCharged * #it{p}_{T}^{3} / tsallisCharged(1)");
  gr.SetLineColor(kRed);
  gr.SetLineWidth(2);
  for(int i=0; i<300; ++i) {
    const double x = i/100.;
//     const double y = rho(x);
    const double y = rho(x) * x*x*x / rho(1);
    gr.AddPoint(x, y);
  }

  TCanvas cc("cc", "", 1200, 800);
  cc.SetLogy();
  gr.Draw("AL");
  cc.Print("gr.pdf", "pdf");
}

double rho(const double pT) {
  const double mt = std::sqrt(mass * mass + pT * pT);
  const double n = a + b / sqrtSNN;
  const double t = c + d / sqrtSNN;
  const double p0 = n * t;
  const double result = std::pow((1. + mt / p0), -n);
  return result;
}
