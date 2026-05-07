std::map<std::string_view, int> MapTHnSparseAxesIndices(const THnSparse* histo);
void SetTHnSparseAxisRanges(THnSparse* histo, int axisNum, float lo, float hi);
TH1* CutSubHistogram(const TH1* histoIn, double lo, double hi);

const std::string_view signalTypeAxisTitle = "candidates type";
const std::string_view pTBAxisTitle = "#it{p}_{T}^{B} (GeV/#it{c})";

const double loPt{0.};
const double hiPt{20.};

void pt_NPweight_builder(const std::string& fileGenName, const std::string& fileFonllName, const std::string& fileInOutName, const std::string& suffix="cent") {
  TFile* fileGen = TFile::Open(fileGenName.c_str(), "read");
  TFile* fileFonll = TFile::Open(fileFonllName.c_str(), "read");

  THnSparse* histoGenTHn = fileGen->Get<THnSparse>("hf-task-lc/hnLcVarsGen");
  const std::map<std::string_view, int> axesIndices = MapTHnSparseAxesIndices(histoGenTHn);
  SetTHnSparseAxisRanges(histoGenTHn, axesIndices.at(signalTypeAxisTitle), 2., 3.);
  TH1* histoGen = histoGenTHn->Projection(axesIndices.at(pTBAxisTitle));

  TH1* histoFonll = fileFonll->Get<TH1>(("hPtFONLLB" + suffix).c_str());

  histoGen = CutSubHistogram(histoGen, loPt, hiPt);
  histoFonll = CutSubHistogram(histoFonll, loPt, hiPt);

//   histoGen->Rebin(10);

  const double integralGen = histoGen->Integral();
  const double integralFonll = histoFonll->Integral();

  histoGen->Scale(1./integralGen);
  histoFonll->Scale(1./integralFonll);

  const bool ok = histoFonll->Divide(histoGen);
  if(!ok) throw std::runtime_error("histoFonll->Divide(histoGen) was not ok");

  TFile* fileOut = TFile::Open(fileInOutName.c_str(), "update");
  histoFonll->Write(("histoNPWeight" + suffix).c_str());
  fileOut->Close();

  fileFonll->Close();
  fileGen->Close();
}


std::map<std::string_view, int> MapTHnSparseAxesIndices(const THnSparse* histo) {
  std::map<std::string_view, int> result;
  const int nDims = histo->GetNdimensions();
  for(int iDim=0; iDim<nDims; ++iDim) {
    result.insert({histo->GetAxis(iDim)->GetTitle(), iDim});
  }
  return result;
}

void SetTHnSparseAxisRanges(THnSparse* histo, int axisNum, float lo, float hi) {
  constexpr double tolerance = 1e-6;

  if(std::fabs(lo+999)<tolerance && std::fabs(hi+999)<tolerance) {
    histo->GetAxis(axisNum)->SetRange();
    return;
  }

  if(lo >= hi) throw std::runtime_error("SetTHnSparseAxisRanges(): lo >= hi");

  const TAxis* axis = histo->GetAxis(axisNum);
  int binLo{-999}, binHi{-999};
  for(int iBin=1, nBins=axis->GetNbins(); iBin<=nBins; ++iBin) {
    const float binLowEdge = axis->GetBinLowEdge(iBin);
    const float binUpEdge = axis->GetBinUpEdge(iBin);
    if(std::fabs(binLowEdge - lo)<tolerance) binLo = iBin;
    if(std::fabs(binUpEdge - hi)<tolerance) binHi = iBin;
    if(binLo != -999 && binHi != -999) break;
  }
  if(binLo == -999 || binHi == -999) throw std::runtime_error("SetTHnSparseAxisRanges(): binLo == -999 || binHi == -999");
  histo->GetAxis(axisNum)->SetRange(binLo, binHi);
}

TH1* CutSubHistogram(const TH1* histoIn, double lo, double hi) {
  if(lo >= hi) throw std::runtime_error("HelperMath::CutSubHistogram(): lo >= hi");

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
  if(binLoIn == -999 || !isEndReached) throw std::runtime_error("HelperMath::CutSubHistogram(): either lo or hi does not match any of histoIn bin edges");

  TH1* histoOut = new TH1D("", "", binEdges.size()-1, binEdges.data());
  histoOut->SetDirectory(nullptr);
  if(histoIn->GetSumw2N() > 0) histoOut->Sumw2();
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

  return histoOut;
}
