TH1* CutSubHistogram(const TH1* histoIn, double lo, double hi);
TF1* FitLifetimeHisto(TH1* histo, const std::string& option);
std::pair<float, float> EstimateExpoParameters(TH1* h, float lo, float hi);
TH1* EvaluateEfficiency(const TH1* histoNum, const TH1* histoDen);

void smearToyMc() {
  gStyle->SetHistLineWidth(2);
  TFile* fileReso = TFile::Open("grRes.vsMc.root", "read");
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

  const int nBins{400};
  const double lo{0.};
  const double hi{2.};
  const double binWidth = (hi - lo) / nBins;
  const std::vector<double> edges{0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8};

  const double tau{0.2};

  const int nFills{100000000};

  TH1* hExpo = new TH1D("hExpo", "", nBins, lo, hi);
  TH1* hSmeared = new TH1D("hSmeared", "", nBins, lo, hi);

  for(int iFill=0; iFill<nFills; ++iFill) {
    if(iFill%(nFills/20) == 0) std::cout << "iFill = " << iFill << "\n";
    hExpo->Fill(gRandom->Exp(tau));
    const double smearCentralValue = gRandom->Exp(tau);
    if(gRandom->Uniform(1) > hEffSim->GetBinContent(hEffSim->FindBin(smearCentralValue))) continue;
    double mean{};
    double sigma{};
    if(smearCentralValue < grResoSigma->GetX()[0]) {
      mean = grResoMean->GetY()[0];
      sigma = grResoSigma->GetY()[0];
    } else if(smearCentralValue > grResoSigma->GetX()[grResoSigma->GetN()-1]) {
      mean = grResoMean->GetY()[grResoMean->GetN()-1];
      sigma = grResoSigma->GetY()[grResoSigma->GetN()-1];
    } else {
      mean = grResoMean->Eval(smearCentralValue);
      sigma = grResoSigma->Eval(smearCentralValue);
    }
    const double smearShiftValue = gRandom->Gaus(mean, sigma);
    hSmeared->Fill(smearCentralValue + smearShiftValue);
  }

  hExpo->SetLineColor(kBlue);
  hExpo->Draw();

  hSmeared->SetLineColor(kRed);
//   hSmeared->Draw("HIST same");

//   TH1* hEffCand = EvaluateEfficiency(hYieldCand, hYieldGen);
//   hSmeared->Sumw2();
//   hSmeared->Divide(hEffCand);
//   hSmeared->Draw("same");


  auto CutSubHistogramL = [&] (TH1* histo) { return CutSubHistogram(histo, edges.front(), edges.back()); };
  TH1* hExpoCut = CutSubHistogramL(hExpo);
  TH1* hSmearedCut = CutSubHistogramL(hSmeared);
  TH1* hYieldCandCut = CutSubHistogramL(hYieldCand);
  TH1* hYieldGenCut = CutSubHistogramL(hYieldGen);

  auto RebinL = [&] (TH1*& histo) {
    histo = dynamic_cast<TH1*>(histo->Rebin(edges.size() - 1,histo->GetName(),edges.data()));
    histo->Scale(1., "width");
  };

  RebinL(hExpoCut);
  RebinL(hSmearedCut);
  RebinL(hYieldCandCut);
  RebinL(hYieldGenCut);

  TH1* hEffCandCut = EvaluateEfficiency(hYieldCandCut, hYieldGenCut);

  hSmearedCut->Sumw2();
  hSmearedCut->Divide(hEffCandCut);

  const std::string option{"I"};
  TF1* fitExpo = FitLifetimeHisto(hExpoCut, option.c_str());
  TF1* fitSmeared = FitLifetimeHisto(hSmearedCut, option.c_str());
  fitExpo->SetLineColor(kBlue);
  fitSmeared->SetLineColor(kRed);

  hExpoCut->SetLineColor(kBlue);
  hExpoCut->Draw();
  hSmearedCut->SetLineColor(kRed);
  hSmearedCut->Draw("same");
  fitExpo->Draw("same");
  fitSmeared->Draw("same");

//   fileReso->Close();
}

TH1* CutSubHistogram(const TH1* histoIn, double lo, double hi) {
  if(lo >= hi) throw std::runtime_error("CutSubHistogram(): lo >= hi");

  const double tolerance = 1e-6;
  int binLoIn{-999};
  bool isEndReached{false};
  std::vector<double> binEdges;
  for(int iBin=1, nBins=histoIn->GetNbinsX(); iBin<=nBins+1; ++iBin) {
    const double binLowEdge = histoIn->GetBinLowEdge(iBin);
    if(std::fabs(binLowEdge - lo) < tolerance) binLoIn = iBin;
    if(binLoIn != -999) binEdges.emplace_back(binLowEdge);
    if(std::fabs(binLowEdge - hi) < tolerance) {
      isEndReached = true;
      break;
    }
  } // histoIn bins
  if(binLoIn == -999 || !isEndReached) throw std::runtime_error("CutSubHistogram(): either lo or hi does not match any of histoIn bin edges");

  TH1* histoOut = new TH1D("", "", binEdges.size()-1, binEdges.data());
  histoOut->SetDirectory(nullptr);
  const bool isSumw2 = histoIn->GetSumw2N() > 0;
  histoOut->GetXaxis()->SetTitle(histoIn->GetXaxis()->GetTitle());
  histoOut->GetYaxis()->SetTitle(histoIn->GetYaxis()->GetTitle());
  histoOut->SetName(histoIn->GetName());
  histoOut->SetTitle(histoIn->GetTitle());
  for(int iBin=1, nBins=binEdges.size()-1; iBin<=nBins; ++iBin) {
    const double value = histoIn->GetBinContent(binLoIn-1 + iBin);
    const double error = histoIn->GetBinError(binLoIn-1 + iBin);
    histoOut->SetBinContent(iBin, value);
    histoOut->SetBinError(iBin, error);
  }
  histoOut->Sumw2(isSumw2);
  return histoOut;
}

TF1* FitLifetimeHisto(TH1* histo, const std::string& option) {
  const double lo = histo->GetBinLowEdge(1) + 1e-3;
  const double hi = histo->GetBinLowEdge(histo->GetNbinsX()+1) - 1e-3;
  const auto parEst = EstimateExpoParameters(histo, lo, hi);
  TF1* fitFunc = new TF1("fitFunc", "[0]*TMath::Exp(-x/[1])", lo, hi);
  fitFunc->SetParameters(parEst.first, parEst.second);
  histo->Fit(fitFunc, ("0"+option).c_str(), "", lo, hi);
  fitFunc->SetLineColor(histo->GetLineColor());

  return fitFunc;
}

std::pair<float, float> EstimateExpoParameters(TH1* h, float lo, float hi) {
  const int ilo = h->FindBin(lo);
  const int ihi = h->FindBin(hi);
  const float flo = h->GetBinContent(ilo)/* * h->GetBinWidth(ilo)*/;
  const float fhi = h->GetBinContent(ihi)/* * h->GetBinWidth(ihi)*/;
  const float tau = (hi-lo)/std::log(flo/fhi);
  const float A = flo / std::exp(-lo/tau);
  return std::make_pair(A, tau);
}

TH1* EvaluateEfficiency(const TH1* histoNum, const TH1* histoDen) {
  TH1* histoEff = dynamic_cast<TH1*>(histoNum->Clone());
  histoEff->Sumw2();
  histoEff->Divide(histoDen);

  return histoEff;
}
