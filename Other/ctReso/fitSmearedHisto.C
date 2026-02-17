std::pair<double, double> EstimateExpoParameters(TH1* h, double lo, double hi);

class ExpoSmearFunction {
public:
  void Init() {
    func_orig_ = new TF1("fOrig", "[0]*TMath::Exp(-x/[1])", 0., 2.);
    func_smear_ = new TF1("fSmear", "gaus(0)", -smear_range_, smear_range_);
  }

  void SetSigmaContainer(TGraph* graph) { graph_sigma_container_ = graph; }

  double operator()(double* xx, double* par) {
    double result{};
    const double x = xx[0];
    double sigma{};
    if(x<graph_sigma_container_->GetX()[0]) {
      sigma = graph_sigma_container_->GetY()[0];
    } else if (x>graph_sigma_container_->GetX()[graph_sigma_container_->GetN()-1]) {
      sigma = graph_sigma_container_->GetY()[graph_sigma_container_->GetN()-1];
    } else {
      sigma = graph_sigma_container_->Eval(x);
    }
    func_smear_->SetParameters(1/std::sqrt(2*3.1415)/sigma, 0, sigma);
    func_orig_->SetParameters(par[0], par[1]);
    const double step = 2 * smear_range_ / (numer_of_smearing_points_ - 1);
    for(int iPoint=0; iPoint<numer_of_smearing_points_; ++iPoint) {
      const double diffFromCenter = -smear_range_ + iPoint*step;
      const double funcValue = func_orig_->Eval(x + diffFromCenter);
      const double smearFactor = func_smear_->Eval(diffFromCenter);
      result += funcValue * smearFactor;
    }
    result *= step;
    return result;
  }

private:
  TF1* func_orig_{nullptr};
  TGraph* graph_sigma_container_{nullptr};
  TF1* func_smear_{nullptr};
  double smear_range_{0.2};
  int numer_of_smearing_points_{1000};
};

void fitSmearedHisto() {
  gStyle->SetHistLineWidth(2);

  TFile* fileInHisto = TFile::Open("hSmeared.fine.root", "read");
  if(fileInHisto == nullptr) throw std::runtime_error("fileInHisto == nullptr");
  TH1* histoIn = fileInHisto->Get<TH1>("hSmeared");
  if(histoIn == nullptr) throw std::runtime_error("histoIn == nullptr");

  TFile* fileReso = TFile::Open("grRes.root", "read");
  if(fileReso == nullptr) throw std::runtime_error("fileReso == nullptr");
  TGraph* grReso = fileReso->Get<TGraph>("grRes");
  if(grReso == nullptr) throw std::runtime_error("grReso == nullptr");

  const double lo = histoIn->GetBinLowEdge(1) + 1e-3;
  const double hi = histoIn->GetBinLowEdge(histoIn->GetNbinsX()+1) - 1e-3;
  const auto parEst = EstimateExpoParameters(histoIn, lo, hi);
//   TF1* fitFunc = new TF1("fitFunc", "[0]*TMath::Exp(-x/[1])", lo, hi);
//   fitFunc->SetParameters(parEst.first, parEst.second);

  ExpoSmearFunction fitFunctorSmeared;
  fitFunctorSmeared.SetSigmaContainer(grReso);
  fitFunctorSmeared.Init();

  TF1* fitFuncSmeared = new TF1("fitFuncSmeared", fitFunctorSmeared, 0., 2., 2);
  fitFuncSmeared->SetParameters(parEst.first, parEst.second);

  histoIn->Fit(fitFuncSmeared, "I", "", lo, hi);

  fitFuncSmeared->SaveAs("f1.root");

//   histoIn->Fit(fitFunc, "I", "", lo, hi); // Usual fitting

}

std::pair<double, double> EstimateExpoParameters(TH1* h, double lo, double hi) {
  const int ilo = h->FindBin(lo);
  const int ihi = h->FindBin(hi);
  const double flo = h->GetBinContent(ilo)/* * h->GetBinWidth(ilo)*/;
  const double fhi = h->GetBinContent(ihi)/* * h->GetBinWidth(ihi)*/;
  const double tau = (hi-lo)/std::log(flo/fhi);
  const double A = flo / std::exp(-lo/tau);
  return std::make_pair(A, tau);
}
