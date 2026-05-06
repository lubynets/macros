template<typename T>
void checkNullptr(T* t) {
  if(t != nullptr) std::cout << "checkNullptr success\n";
  else throw;
}

std::map<std::string_view, int> MapTHnSparseAxesIndices(const THnSparse* histo);
void SetTHnSparseAxisRanges(THnSparse* histo, int axisNum, float lo, float hi);
TH1* CutSubHistogram(const TH1* histoIn, double lo, double hi);

const std::string_view signalTypeAxisTitle = "candidates type";
const std::string_view pTBAxisTitle = "#it{p}_{T}^{B} (GeV/#it{c})";

const double loPt{0.};
const double hiPt{20.};

void plotNPPtshapes() {
  TFile* fileFonll = TFile::Open("/home/oleksii/alidir/working/tsallis/nonprompt/PtWeigths_NonPromptLc_LHC20l3a_Lb.root");
  TFile* fileMcBWrong = TFile::Open("/home/oleksii/alidir/working/tsallis/HF_LHC24h1b_All/ptGenNonPrompt.HF_LHC24h1b_All.595984.root");
  TFile* fileMcBCorrect = TFile::Open("/home/oleksii/alidir/working/cutVar/mcClosure/input/AnalysisResults.HL.mc.HF_LHC24h1b_All.667716.root");
  TFile* fileMcC = TFile::Open("/home/oleksii/alidir/working/cutVar/mcClosure/input/AnalysisResults.HL.mc.HF_LHC24h1c_All.669133.root");

  checkNullptr(fileFonll);
  checkNullptr(fileMcBWrong);
  checkNullptr(fileMcBCorrect);
  checkNullptr(fileMcC);

  TH1* hFonll = fileFonll->Get<TH1>("hPtFONLLBcent");
  TH1* hWrong = fileMcBWrong->Get<TH1>("hPtLambdaB");
  THnSparse* hThnCorrect = fileMcBCorrect->Get<THnSparse>("hf-task-lc/hnLcVarsGen");
  THnSparse* hThnNew = fileMcC->Get<THnSparse>("hf-task-lc/hnLcVarsGen");

  checkNullptr(hFonll);
  checkNullptr(hWrong);
  checkNullptr(hThnCorrect);
  checkNullptr(hThnNew);

  const std::map<std::string_view, int> axesIndices = MapTHnSparseAxesIndices(hThnCorrect);
  SetTHnSparseAxisRanges(hThnCorrect, axesIndices.at(signalTypeAxisTitle), 2., 3.);
  SetTHnSparseAxisRanges(hThnNew, axesIndices.at(signalTypeAxisTitle), 2., 3.);

  TH1* hCorrect = hThnCorrect->Projection(axesIndices.at(pTBAxisTitle));
  TH1* hNew = hThnNew->Projection(axesIndices.at(pTBAxisTitle));

  checkNullptr(hCorrect);
  checkNullptr(hNew);

  hWrong->Rebin(10);
  for(auto& histo : {&hFonll, &hWrong, &hCorrect, &hNew}) {
    (*histo) = CutSubHistogram(*histo, loPt, hiPt);
    (*histo)->Scale(1./(*histo)->Integral());
    (*histo)->SetLineWidth(3);
  }

  hFonll->SetLineColor(kRed);
  hWrong->SetLineColor(kGreen+2);
  hCorrect->SetLineColor(kMagenta);
  hNew->SetLineColor(kBlue);

  TCanvas* cc = new TCanvas("cc", "", 1200, 800);
  hFonll->Draw();
  hWrong->Draw("same");
  hCorrect->Draw("same");
  hNew->Draw("same");
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
