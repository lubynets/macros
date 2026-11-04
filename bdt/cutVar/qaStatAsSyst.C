template<std::size_t N>
double standard_deviation(const std::array<double, N>& a);

void qaStatAsSyst() {
  const std::array<std::string, 3> sigmas{"Min", "Max", "Cent"};

  std::array<std::vector<double>, sigmas.size()*sigmas.size()> values{};

  std::vector<std::array<double, sigmas.size()*sigmas.size()>> Values{};

  int iSigma{0};
  TH1* hSyst{nullptr};
  for(const auto& sigmaP : sigmas) {
    for(const auto& sigmaNP : sigmas) {
      TFile* fileIn = TFile::Open(("/home/oleksii/alidir/working/cutVar/testZeroEff/" + sigmaP + "/" + sigmaNP + "/CutVarLc.merged.root").c_str(), "");
      TH1* histoIn = fileIn->Get<TH1>("hCorrYieldsPrompt");
      if(iSigma == 0) {
        hSyst = dynamic_cast<TH1*>(histoIn->Clone());
        hSyst->SetDirectory(nullptr);
        hSyst->Reset();

        Values.resize(histoIn->GetNbinsX());
      }
      for(int iBin=1, nBins=histoIn->GetNbinsX(); iBin<=nBins; ++iBin) {
        values.at(iSigma).push_back(histoIn->GetBinContent(iBin));
        Values.at(iBin-1).at(iSigma) = histoIn->GetBinContent(iBin);
      }
      ++iSigma;
      fileIn->Close();
    }
  }

  for(int iBin=1, nBins=hSyst->GetNbinsX(); iBin<=nBins; ++iBin) {
    hSyst->SetBinContent(iBin, standard_deviation(Values.at(iBin-1)));
  }

  hSyst->SaveAs("hSyst.root");
}

template<std::size_t N>
double standard_deviation(const std::array<double, N>& a) {
    static_assert(N > 1, "Need at least two measurements");

    const double mean =
        std::accumulate(a.begin(), a.end(), 0.0) / N;

    const double squared_deviations =
        std::accumulate(a.begin(), a.end(), 0.0,
            [mean](double sum, double x) {
                const double d = x - mean;
                return sum + d * d;
            });

    return std::sqrt(squared_deviations / (N - 1));
}
