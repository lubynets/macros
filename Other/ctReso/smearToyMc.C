TH1* CutSubHistogram(const TH1* histoIn, double lo, double hi);
TF1* FitLifetimeHisto(TH1* histo, const std::string& option);
std::pair<float, float> EstimateExpoParameters(TH1* h, float lo, float hi);

void smearToyMc() {
  gStyle->SetHistLineWidth(2);
  TFile* fileReso = TFile::Open("grRes.root", "read");
  if(fileReso == nullptr) throw std::runtime_error("fileReso == nullptr");
  TGraph* grReso = fileReso->Get<TGraph>("grRes");
  if(grReso == nullptr) throw std::runtime_error("grReso == nullptr");

  const int nBins{400};
  const double lo{0.};
  const double hi{2.};
  const double binWidth = (hi - lo) / nBins;
  const std::vector<double> edges{0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8};

  const double tau = 0.203729;

  const int nFills{100000000};

  TH1* hExpo = new TH1D("hExpo", "HEXPO", nBins, lo, hi);
  TH1* hSmeared = new TH1D("hSmeared", "HSMEARED", nBins, lo, hi);

  for(int iFill=0; iFill<nFills; ++iFill) {
    hExpo->Fill(gRandom->Exp(tau));
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

  hExpo->SetLineColor(kBlue);
//   hExpo->Draw();

  const double smearRange{1.};

//   const double smearSigma{0.2};
//   TF1* fGaus = new TF1("fGaus", "gaus(0)", -smearRange, smearRange);
//   fGaus->SetParameters(1/std::sqrt(2*3.1415)/smearSigma, 0, smearSigma);
// //   fGaus->Draw();

//   TF1* fGaus = new TF1("fGaus", [&] (Double_t *x, Double_t *p) {
//     const double sigma = grReso->Eval(x[0]);
//     return TMath::Gaus(x[0], 0, sigma, true);
//   }, -smearRange, smearRange, 0);
//
//   for(int iBinTarget=1; iBinTarget<=nBins; ++iBinTarget) {
//     for(int iBinSource=1; iBinSource<=nBins; ++iBinSource) {
//       const double centerTarget = hExpo->GetBinCenter(iBinTarget);
//       const double centerSource = hExpo->GetBinCenter(iBinSource);
//       const double diff = centerSource - centerTarget;
//       if(std::fabs(diff) > smearRange) continue;
//       const double weight = hExpo->GetBinContent(iBinSource) * fGaus->Eval(diff);
//       hSmeared->Fill(centerTarget, weight);
//     }
//   }
//   hSmeared->Scale(binWidth);
//
//   for(int iBinTarget=1; iBinTarget<=nBins; ++iBinTarget) {
//     const double value = hSmeared->GetBinContent(iBinTarget);
//     hSmeared->SetBinError(iBinTarget, std::sqrt(value));
//   }
//   hSmeared->Sumw2(false);




  hSmeared->SetLineColor(kRed);
//   hSmeared->Draw("HIST same");

//   TH1* hRatio = dynamic_cast<TH1*>(hSmeared->Clone());
//   hRatio->Sumw2();
//   hRatio->Divide(hExpo);
//   hRatio->Draw();

//
  TH1* hExpoCut = CutSubHistogram(hExpo, edges.front(), edges.back());
  TH1* hSmearedCut = CutSubHistogram(hSmeared, edges.front(), edges.back());

//   hExpoCut = dynamic_cast<TH1*>(hExpoCut->Rebin(edges.size() - 1,hExpoCut->GetName(),edges.data()));
//   hSmearedCut = dynamic_cast<TH1*>(hSmearedCut->Rebin(edges.size() - 1,hSmearedCut->GetName(),edges.data()));
//
//   hExpoCut->Scale(1., "width");
//   hSmearedCut->Scale(1., "width");

  hSmearedCut->SaveAs("hSmeared.root");

//   TF1* fitExpo = FitLifetimeHisto(hExpoCut, "I");
//   TF1* fitSmeared = FitLifetimeHisto(hSmearedCut, "I");
//   fitExpo->SetLineColor(kBlue);
//   fitSmeared->SetLineColor(kRed);
//
//   hExpoCut->SetLineColor(kBlue);
//   hExpoCut->Draw();
//   hSmearedCut->SetLineColor(kRed);
//   hSmearedCut->Draw("same");
//   fitExpo->Draw("same");
//   fitSmeared->Draw("same");

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
