constexpr double PI{3.141592654};

TH1* EvaluateEfficiency(const TH1* histoNum, const TH1* histoDen);

class ExpoSmearFunction {
public:
  void Init() {
    func_orig_ = new TF1("fOrig", "x < 0 ? 0. : [0]*TMath::Exp(-x/[1])", 0., 2.);
    func_smear_ = new TF1("fSmear", "gaus(0)", -smear_range_, smear_range_);
    step_ = 2 * smear_range_ / (numer_of_smearing_points_ - 1);
    is_init_ = true;
  }

  void SetMeanContainer(TGraph* graph) { graph_mean_container_ = graph; }
  void SetSigmaContainer(TGraph* graph) { graph_sigma_container_ = graph; }
  void SetEfficiencyGenHistogram(TH1* histo) { efficiency_gen_based_ = histo; }
  void SetEfficiencyRecHistogram(TH1* histo) { efficiency_rec_based_ = histo; }

  double operator()(double* xx, double* par) {
    if(!is_init_) Init();
    double result{};
    const double x = xx[0];

    func_orig_->SetParameters(par[0], par[1]);

    auto ReadValueFromGraph = [](TGraph* graph, double arg) {
      if(arg < graph->GetX()[0]) {
        return  graph->GetY()[0];
      } else if(arg > graph->GetX()[graph->GetN()-1]) {
        return  graph->GetY()[graph->GetN()-1];
      } else {
        return  graph->Eval(arg);
      }
    };

    for(int iPoint=0; iPoint<numer_of_smearing_points_; ++iPoint) {
      const double shiftFromPointToBeFilled = -smear_range_ + iPoint*step_;
      const double pointToBeTakenFrom = x + shiftFromPointToBeFilled;

      const double mean = ReadValueFromGraph(graph_mean_container_, pointToBeTakenFrom);
      const double sigma = ReadValueFromGraph(graph_sigma_container_, pointToBeTakenFrom);
      const double effGen = efficiency_gen_based_->Interpolate(pointToBeTakenFrom);
      const double effRec = efficiency_rec_based_->Interpolate(x);

      const double funcValue = func_orig_->Eval(pointToBeTakenFrom);

      func_smear_->SetParameters(1/std::sqrt(2*PI)/sigma, 0, sigma);

      const double smearFactor = func_smear_->Eval(shiftFromPointToBeFilled + mean);
      result += funcValue * smearFactor * effGen / effRec;
    }
    result *= step_;
    return result;
  }

private:
  TF1* func_orig_{nullptr};
  TGraph* graph_mean_container_{nullptr};
  TGraph* graph_sigma_container_{nullptr};
  TF1* func_smear_{nullptr};
  TH1* efficiency_gen_based_{nullptr};
  TH1* efficiency_rec_based_{nullptr};
  double smear_range_{1.};
  double step_{};
  int numer_of_smearing_points_{1000};
  bool is_init_{false};
};

void th1tf1SmearConsistency() {
  gStyle->SetHistLineWidth(2);
  TFile* fileReso = TFile::Open("grRes/grRes.prompt.unsigned.root", "read");
  if(fileReso == nullptr) throw std::runtime_error("fileReso == nullptr");
  TGraph* grResoMean = fileReso->Get<TGraph>("mean");
  if(grResoMean == nullptr) throw std::runtime_error("grResoMean == nullptr");
  TGraph* grResoSigma = fileReso->Get<TGraph>("sigma");
  if(grResoSigma == nullptr) throw std::runtime_error("grResoSigma == nullptr");

  TFile* fileYield = TFile::Open("ct_yield_qa.HF_LHC24h1b_All.595984.root", "read");
  if(fileYield == nullptr) throw std::runtime_error("fileYield == nullptr");
  TH1* hYieldGen = fileYield->Get<TH1>("prompt/hGen");
  TH1* hYieldCand = fileYield->Get<TH1>("prompt/hCand");
  TH1* hYieldSim = fileYield->Get<TH1>("prompt/hSim");
  if(hYieldGen == nullptr || hYieldCand == nullptr || hYieldSim == nullptr) throw std::runtime_error("hYieldGen == nullptr || hYieldCand == nullptr || hYieldSim == nullptr");

  TH1* hEffSim = EvaluateEfficiency(hYieldSim, hYieldGen);
  TH1* hEffCand = EvaluateEfficiency(hYieldCand, hYieldGen);

  const int nBins{400};
  const double lo{0.};
  const double hi{2.};
  const double binWidth = (hi - lo) / nBins;

  const double tau{0.2};

  const int nFills{10000000};
  TH1* hSmeared = new TH1D("hSmeared", "", nBins, lo, hi);
  for(int iFill=0; iFill<nFills; ++iFill) {
    if(iFill%(nFills/20) == 0) std::cout << "iFill = " << iFill << "\n";
    const double smearCentralValue = gRandom->Exp(tau);
    if(gRandom->Uniform(1) > hEffSim->GetBinContent(hEffSim->FindBin(smearCentralValue))) continue;

    auto ReadValueFromGraph = [&](TGraph* graph) {
      if(smearCentralValue < graph->GetX()[0]) {
        return  graph->GetY()[0];
      } else if(smearCentralValue > graph->GetX()[graph->GetN()-1]) {
        return  graph->GetY()[graph->GetN()-1];
      } else {
        return  graph->Eval(smearCentralValue);
      }
    };
    const double mean = ReadValueFromGraph(grResoMean);
    const double sigma = ReadValueFromGraph(grResoSigma);

    const double smearShiftValue = gRandom->Gaus(0, sigma);
    hSmeared->Fill(smearCentralValue + smearShiftValue);
  }

  hSmeared->Sumw2();
  hSmeared->Divide(hEffCand);

  hSmeared->SetLineColor(kBlue);

  //================================================================================================

  ExpoSmearFunction fitFunctorSmeared;
  fitFunctorSmeared.SetMeanContainer(grResoMean);
  fitFunctorSmeared.SetSigmaContainer(grResoSigma);
  fitFunctorSmeared.SetEfficiencyGenHistogram(hEffSim);
  fitFunctorSmeared.SetEfficiencyRecHistogram(hEffCand);
  fitFunctorSmeared.Init();

  TF1* fitFuncSmeared = new TF1("fitFuncSmeared", fitFunctorSmeared, 0.0, 2.0, 2);
  fitFuncSmeared->SetNpx(1000);

  const double A = nFills * binWidth / tau / (std::exp(-lo/tau) - std::exp(-hi/tau));
  fitFuncSmeared->SetParameters(A, tau);

//   hSmeared->Draw();
  fitFuncSmeared->Draw();

//   hSmeared->Divide(fitFuncSmeared);
//   hSmeared->Draw();

}

TH1* EvaluateEfficiency(const TH1* histoNum, const TH1* histoDen) {
  TH1* histoEff = dynamic_cast<TH1*>(histoNum->Clone());
  histoEff->Sumw2();
  histoEff->Divide(histoDen);

  return histoEff;
}
