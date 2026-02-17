class ExpoSmearFunction {
public:
  void Init() {
    func_orig_ = new TF1("fOrig", "x < 0 ? 0. : [0]*TMath::Exp(-x/[1])", 0., 2.);
    func_smear_ = new TF1("fSmear", "gaus(0)", -smear_range_, smear_range_);
    step_ = 2 * smear_range_ / (numer_of_smearing_points_ - 1);
    is_init_ = true;
  }

  void SetSigmaContainer(TGraph* graph) { graph_sigma_container_ = graph; }

  double operator()(double* xx, double* par) {
    if(!is_init_) Init();
    double result{};
    const double x = xx[0];
    double sigma{};
    if(x < graph_sigma_container_->GetX()[0]) {
      sigma = graph_sigma_container_->GetY()[0];
    } else if (x > graph_sigma_container_->GetX()[graph_sigma_container_->GetN()-1]) {
      sigma = graph_sigma_container_->GetY()[graph_sigma_container_->GetN()-1];
    } else {
      sigma = graph_sigma_container_->Eval(x);
    }
    func_smear_->SetParameters(1/std::sqrt(2*3.1415)/sigma, 0, sigma);
    func_orig_->SetParameters(par[0], par[1]);
    for(int iPoint=0; iPoint<numer_of_smearing_points_; ++iPoint) {
      const double diffFromCenter = -smear_range_ + iPoint*step_;
      const double funcValue = func_orig_->Eval(x + diffFromCenter);
      const double smearFactor = func_smear_->Eval(diffFromCenter);
      result += funcValue * smearFactor;
    }
    result *= step_;
    return result;
  }

private:
  TF1* func_orig_{nullptr};
  TGraph* graph_sigma_container_{nullptr};
  TF1* func_smear_{nullptr};
  double smear_range_{1.};
  double step_{};
  int numer_of_smearing_points_{1000};
  bool is_init_{false};
};

void th1tf1SmearConsistency() {
  gStyle->SetHistLineWidth(2);
  TFile* fileReso = TFile::Open("grRes.root", "read");
  if(fileReso == nullptr) throw std::runtime_error("fileReso == nullptr");
  TGraph* grReso = fileReso->Get<TGraph>("grRes");
  if(grReso == nullptr) throw std::runtime_error("grReso == nullptr");

  const int nBins{400};
  const double lo{0.};
  const double hi{2.};
  const double binWidth = (hi - lo) / nBins;

  const double tau{0.2};

  const int nFills{100000000};
  TH1* hSmeared = new TH1D("hSmeared", "", nBins, lo, hi);
  for(int iFill=0; iFill<nFills; ++iFill) {
    const double smearCentralValue = gRandom->Exp(tau);
    double sigma{};
    if(smearCentralValue < grReso->GetX()[0]) {
      sigma = grReso->GetY()[0];
    } else if(smearCentralValue > grReso->GetX()[grReso->GetN()-1]) {
      sigma = grReso->GetY()[grReso->GetN()-1];
    } else {
      sigma = grReso->Eval(smearCentralValue);
    }
    const double smearShiftValue = gRandom->Gaus(0, sigma);
    hSmeared->Fill(smearCentralValue + smearShiftValue);
  }

  hSmeared->SetLineColor(kBlue);

  //================================================================================================

  ExpoSmearFunction fitFunctorSmeared;
  fitFunctorSmeared.SetSigmaContainer(grReso);
  fitFunctorSmeared.Init();

  TF1* fitFuncSmeared = new TF1("fitFuncSmeared", fitFunctorSmeared, 0., 2., 2);
  fitFuncSmeared->SetNpx(1000);

  const double A = nFills * binWidth / tau / (std::exp(-lo/tau) - std::exp(-hi/tau));
  fitFuncSmeared->SetParameters(A, tau);

//   hSmeared->Draw();
//   fitFuncSmeared->Draw("same");

  hSmeared->Divide(fitFuncSmeared);
  hSmeared->Draw();

}
