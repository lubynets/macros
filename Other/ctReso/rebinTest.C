std::pair<double, double> EstimateExpoParameters(TH1* h, double lo, double hi);
double mean(const std::vector<double>& v);
double sampleStdDev(const std::vector<double>& v);
double standardErrorOfMean(const std::vector<double>& v);

void rebinTest() {
  const int nTrials{10};
  const int nFills{100000};
  const double tau{10.};
  const int nBins{60};
  const double lo{0.};
  const double hi{60.};

  auto doFit = [&](const int rebinFactor = 1) {

    std::vector<double> errors{}, residuals{}, pulls{};
    TH1* histoPull = new TH1F("histoPull", "", 100, -5., 5.);

    for(int iTrial=0; iTrial<nTrials; ++iTrial) {
      TH1* histo = new TH1F("histo", "", nBins, lo, hi);
      for(int iFill=0; iFill<nFills; ++iFill) {
        histo->Fill(gRandom->Exp(tau));
      }
      if(rebinFactor != 1) histo->Rebin(rebinFactor);
      TF1* fitFunc = new TF1("fitFunc", "[0]*TMath::Exp(-x/[1])", lo, hi);
      const auto parEst = EstimateExpoParameters(histo, lo, hi);
      fitFunc->SetParameters(parEst.first, parEst.second);
      histo->Fit(fitFunc, "0ILQ", "", lo, hi);
      const double value = fitFunc->GetParameter(1);
      const double error = fitFunc->GetParError(1);
      errors.push_back(error);
      residuals.push_back(value - tau);
      pulls.push_back(residuals.back() / error);
      histoPull->Fill(pulls.back());
      delete histo;
      delete fitFunc;
    }

    const double meanOfError = mean(errors);
    const double meanOfErrorError = standardErrorOfMean(errors);

  //   histoPull->Draw();

    std::cout << "\n\nmeanOfError = " << meanOfError << "\n";
    std::cout << "pull.mean = " << histoPull->GetMean() << " +- " << histoPull->GetMeanError() << "\n";
    std::cout << "pull.stddev = " << histoPull->GetStdDev() << " +- " << histoPull->GetStdDevError() << "\n";

    delete histoPull;
    return std::make_pair(meanOfError, meanOfErrorError);
  };

  const std::vector<int> rebinFactors{1, 2, 3, 4, 5, 6, 10, 12, 15, 20};

  TGraphErrors* gr = new TGraphErrors();
  gr->SetMarkerStyle(kFullSquare);
  gr->SetLineWidth(2);
  gr->GetXaxis()->SetTitle("# of bins");
  gr->GetYaxis()->SetTitle("average error of tau");

  for(const auto& rf : rebinFactors) {
    const auto [meanOfErr, meanOfErrErr] = doFit(rf);
    gr->AddPoint(nBins / rf, meanOfErr);
    gr->SetPointError(gr->GetN()-1, 0., meanOfErrErr);
  }

  gr->Draw("AP");
}

std::pair<double, double> EstimateExpoParameters(TH1* h, double lo, double hi) {
  const int ilo = h->FindBin(lo + 0.001);
  const int ihi = h->FindBin(hi - 0.001);
  const double flo = h->GetBinContent(ilo)/* * h->GetBinWidth(ilo)*/;
  const double fhi = h->GetBinContent(ihi)/* * h->GetBinWidth(ihi)*/;
  const double tau = (hi-lo)/std::log(flo/fhi);
  const double A = flo / std::exp(-lo/tau);
  return std::make_pair(A, tau);
}

double mean(const std::vector<double>& v) {
  return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}

double sampleStdDev(const std::vector<double>& v) {
  double m = mean(v);

  double sqsum = 0.0;
  for (double x : v) {
    double d = x - m;
    sqsum += d * d;
  }

  return std::sqrt(sqsum / (v.size() - 1));
}

double standardErrorOfMean(const std::vector<double>& v) {
  return sampleStdDev(v) / std::sqrt(v.size());
}
