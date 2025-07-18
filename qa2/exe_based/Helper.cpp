//
// Created by oleksii on 11.02.25.
//
#include "Helper.hpp"

#include <TCanvas.h>
#include <TROOT.h>

void Helper::SlightlyShiftXAxis(TGraph* gr, float value) {
  if(value == -1) {
    value = (gr->GetPointX(gr->GetN()-1) - gr->GetPointX(0)) / 50;
  }
  const int nPoints = gr->GetN();
  for(int iPoint=0; iPoint<nPoints; iPoint++) {
    gr->SetPointX(iPoint, gr->GetPointX(iPoint) + value);
  }
}

Helper::HistoQuantities Helper::EvaluateHistoQuantities(const TH1* h) {
  HistoQuantities result;
  const float integral = h->Integral(0, h->GetNbinsX()+1);
  result.nentries_ = integral;
  result.underflow_ = h->GetBinContent(0) / integral;
  result.overflow_ = h->GetBinContent(h->GetNbinsX() + 1) / integral;
  result.mean_ = h->GetMean();
  result.mean_err_ = h->GetMeanError();
  result.stddev_ = h->GetStdDev();
  result.stddev_err_ = h->GetStdDevError();

  return result;
}

TPaveText Helper::ConvertHistoQuantitiesToText(const HistoQuantities& q, float x1, float y1, float x2, float y2) {
  TPaveText text(x1, y1, x2, y2, "brNDC");
  text.SetFillColor(0);
  text.SetTextSize(0.03);
  text.SetTextFont(62);

  text.AddText(("nentries = " + to_string_with_precision(q.nentries_, 0)).c_str());
  text.AddText(("underflow = " + to_string_with_precision(q.underflow_ * 100, 2) + "%").c_str());
  text.AddText(("overflow = " + to_string_with_precision(q.overflow_ * 100, 2) + "%").c_str());
  text.AddText(
          ("#mu = " + to_string_with_precision(q.mean_, 3) + " #pm " + to_string_with_precision(q.mean_err_, 3) +
           " (stat.)").c_str());
  text.AddText(("#sigma = " + to_string_with_precision(q.stddev_, 3) + " #pm " +
                to_string_with_precision(q.stddev_err_, 3) + " (stat.)").c_str());

  return text;
}

std::string Helper::getSubstringBeforeLastSlash(const std::string& input) {
  size_t pos = input.rfind('/');  // Find the last occurrence of the slash

  if (pos == std::string::npos) {
    return "";  // No slash found, return an empty string
  }

  return input.substr(0, pos);  // Return the substring before the last slash
}

std::vector<std::pair<std::string, std::string>> Helper::FindCuts(TFile* fileIn, std::string name_start, bool printCuts) {
  if (name_start.back() != '_') name_start.push_back('_');
  std::vector<std::pair<std::string, std::string>> result;

  const std::string name_start_before_slash = getSubstringBeforeLastSlash(name_start);
  auto lok = fileIn->GetDirectory(name_start_before_slash.c_str())->GetListOfKeys();
  const int nDirs = lok->GetEntries();
  for (int iDir = 0; iDir < nDirs; iDir++) {
    const std::string dirName = lok->At(iDir)->GetName();
    const std::string fullDirName = !name_start_before_slash.empty() ? name_start_before_slash + "/" + dirName : dirName;
    if (fullDirName.substr(0, name_start.size()) != name_start) continue;
    std::pair<std::string, std::string> cutPair;
    bool isFirstCutRead{false};
    bool isGoodCut{true};
    for (int iChar = name_start.size(); iChar < fullDirName.size(); iChar++) {
      char letter = fullDirName.at(iChar);
      if(!(std::isdigit(letter) || letter == '.' || letter == '_')) {
        isGoodCut = false;
        break;
      }
      if (letter != '_') {
        if (!isFirstCutRead) cutPair.first.push_back(letter);
        else cutPair.second.push_back(letter);
      } else {
        isFirstCutRead = true;
      }
    }
    if(isGoodCut) result.emplace_back(cutPair);
  }

  std::sort(result.begin(), result.end(), stofCompare);

  if (result.empty()) {
    throw std::runtime_error("FindCuts(): " + name_start + " cuts are not present");
  }

  if(printCuts) {
    std::cout << "Slice cuts are:\n";
    for(auto& r : result) std::cout << r.first << "\t" << r.second << "\n";
  }

  return result;
}

bool Helper::stofCompare(std::pair<std::string, std::string> a, std::pair<std::string, std::string> b) {
  return atof(a.first.c_str()) < atof(b.first.c_str());
}

void Helper::CustomizeGraphYRange(TGraphMultiErrors* graph, int ne, TF1* f) {
  const int nPoints = graph->GetN();
  float min = 1e9;
  float max = -1e9;

  for (int ie = 0; ie < ne; ie++) {
    for (int iPoint = 0; iPoint < nPoints; iPoint++) {
      const float up = graph->GetPointY(iPoint) + graph->GetErrorY(iPoint, ie);
      const float lo = graph->GetPointY(iPoint) - graph->GetErrorY(iPoint, ie);

      min = std::min(min, lo);
      max = std::max(max, up);
    }
  }

  if (f != nullptr) {
    const float lineLevel = f->GetParameter(0);
    min = std::min(min, lineLevel);
    max = std::max(max, lineLevel);
  }

  const float diff = max - min;
  max += diff / 10;
  min -= diff / 10;

  graph->GetYaxis()->SetRangeUser(min, max);
}

void Helper::SetLineDrawParameters(std::vector<TF1*> fs, int lineWidth, int lineStyle, Color_t lineColor) {
  for (auto& f: fs) {
    f->SetLineWidth(lineWidth);
    f->SetLineStyle(lineStyle);
    f->SetLineColor(lineColor);
  }
}

TF1* Helper::HorizontalLine4Graph(float level, TGraph* graph) {
  float xlo = graph->GetPointX(0);
  float xhi = graph->GetPointX(graph->GetN() - 1);
  const float diff = xhi - xlo;
  xlo = -diff / 10;
  xhi += diff / 10;
  TF1* horizLine = new TF1("horizLine", "[0]", xlo, xhi);
  horizLine->SetParameter(0, level);

  return horizLine;
}

TPaveText* Helper::AddOneLineText(const std::string& text, const std::array<float, 4>& xy, const std::string& option, float size) {
  if(text.empty()) return nullptr;
  TPaveText* textPtr = new TPaveText(xy.at(0), xy.at(1), xy.at(2), xy.at(3), option.c_str());
  textPtr->SetFillColor(0);
  textPtr->SetTextSize(size);
  textPtr->SetTextFont(62);
  textPtr->AddText(text.c_str());
  textPtr->Draw("same");
  return textPtr;
}

std::vector<TPaveText*> Helper::AddMultiLineText(const std::vector<std::string>& texts, const std::array<float, 4>& xy, const std::string& option, float size) {
  const float step = xy.at(3) - xy.at(1);
  std::vector<TPaveText*> ptr;
  for(int iLine=0, nLines=texts.size(); iLine<nLines; iLine++) {
    ptr.emplace_back(AddOneLineText(texts.at(iLine), {xy.at(0), xy.at(1)-iLine*step, xy.at(2), xy.at(3)-iLine*step}, option, size));
  }
  return ptr;
}

bool Helper::string_to_bool(const std::string& str) {
  if(str == "true") return true;
  else if(str == "false") return false;
  else throw std::runtime_error("string_to_bool(): argument must be either true or false");
}

std::pair<double, double> Helper::GetMinMaxBinWithError(const TH1* h) {
  double min = 1e9;
  double max = -1e9;
  for(int iBin=1; iBin<=h->GetNbinsX(); iBin++) {
    min = std::min(min, h->GetBinContent(iBin) - h->GetBinError(iBin));
    max = std::max(max, h->GetBinContent(iBin) + h->GetBinError(iBin));
  }
  return std::make_pair(min, max);
}

TFile* Helper::OpenFileWithNullptrCheck(const std::string& fileName, const std::string& option) {
  TFile* file = TFile::Open(fileName.c_str(), option.c_str());
  if(file == nullptr) {
    throw std::runtime_error("Helper::OpenFileWithNullptrCheck() - file " + fileName + " is missing");
  }
  return file;
}

void Helper::CustomizeHistogramsYRange(const std::vector<TH1*>& histos, bool isLog, double lo, double hi, double part) {
  double max = -1e9;
  double min = 1e9;
  for(auto& histo : histos) {
    auto [hmin, hmax] = GetMinMaxBinWithError(histo);
    max = std::max(max, hmax);
    min = std::min(min, hmin);
  }
  if(isLog) {
    max = std::log10(max);
    min = std::log10(min);
  }
  const double diff = max - min;
  double up = std::min(hi, max + (1-part)/2*diff/part);
  double down = std::max(lo, min - (1-part)/2*diff/part);
  if(isLog) {
    up = std::pow(10, up);
    down = std::pow(10, down);
  }
  for(auto& histo : histos) {
    histo->GetYaxis()->SetRangeUser(down, up);
  }
}

void Helper::PrintInfoOnTF1(const TF1* f) {
  std::cout << "Name = " << f->GetName() << "\n";
  std::cout << "Title = " << f->GetTitle() << "\n";
  const int nPar = f->GetNpar();
  std::cout << "NPar = " << nPar << "\n";
  std::cout << "NFreePar = " << f->GetNumberFreeParameters() << "\n";
  for(int iPar = 0; iPar<nPar; iPar++) {
    double min, max;
    f->GetParLimits(iPar, min, max);
    std::cout << f->GetParName(iPar) << "\t" << f->GetParameter(iPar) << " +- " << f->GetParError(iPar) << "\t(" << min << "; " << max << ")\n";
  }
}

void Helper::CloseCanvasPrinting(const std::vector<std::string>& names) {
  TCanvas emptyCanvas("emptyCanvas", "", 1000, 1000);
  for(auto& name : names) {
    emptyCanvas.Print((name + ".pdf]").c_str(), "pdf");
  }
}

void Helper::CheckHistogramsForXaxisIdentity(const TH1* h1, const TH1* h2) {
  if(h1->GetNbinsX() != h2->GetNbinsX()) {
    throw std::runtime_error("Helper::CheckHistogramsForXaxisIdentity(): nBinsX do not match for " + static_cast<std::string>(h1->GetName()) + " and " + h2->GetName());
  }
  const int nBins = h1->GetNbinsX();
  for(int iBin=1; iBin<=nBins; iBin++) {
    if(std::abs(h1->GetBinCenter(iBin) - h2->GetBinCenter(iBin)) > 1e-6) {
      throw std::runtime_error("Helper::CheckHistogramsForXaxisIdentity(): bins do not coincide for " + static_cast<std::string>(h1->GetName()) + " and " + h2->GetName());
    }
  }
}

std::pair<float, float> Helper::EstimateExpoParameters(TH1* h, float lo, float hi) {
  const int ilo = h->FindBin(lo);
  const int ihi = h->FindBin(hi);
  const float flo = h->GetBinContent(ilo)/* * h->GetBinWidth(ilo)*/;
  const float fhi = h->GetBinContent(ihi)/* * h->GetBinWidth(ihi)*/;
  const float tau = (hi-lo)/std::log(flo/fhi);
  const float A = flo / std::exp(-lo/tau);
  return std::make_pair(A, tau);
}

void Helper::ScaleCanvasVertically(TCanvas* cTo, const TCanvas* cFrom, double scaleFactor) {
  cTo->SetCanvasSize(cFrom->GetWw(), cFrom->GetWh()/scaleFactor);
  cTo->SetBottomMargin(cFrom->GetBottomMargin()*scaleFactor);
  cTo->SetTopMargin(cFrom->GetTopMargin()*scaleFactor);
}

void Helper::LoadMacro(const std::string& macroName) {
  TString currentMacroPath = __FILE__;
  TString directory = currentMacroPath(0, currentMacroPath.Last('/'));
  gROOT->Macro( directory + "/" + macroName );
}

std::pair<double, double> Helper::DetermineWorkingRangesTH1(const TH1* histo, double leftMargin, double rightMargin) {
  double left{-999.};
  double right = {-999.};
  double leftTail = 0.;
  double rightTail = 0.;

  const double integral = histo->Integral(0, histo->GetNbinsX()+1);
  const int nBins = histo->GetNbinsX();

  for(int iBin=0; iBin<=nBins+1; iBin++) {
    leftTail += histo->GetBinContent(iBin);
    if(leftTail > integral*leftMargin) {
      left = iBin == 0 ? -std::numeric_limits<double>::infinity() : histo->GetBinLowEdge(iBin);
      break;
    }
  }

  for(int iBin=nBins+1; iBin>=0; iBin--) {
    rightTail += histo->GetBinContent(iBin);
    if(rightTail > integral*rightMargin) {
      right = iBin == nBins+1 ? std::numeric_limits<double>::infinity() : histo->GetBinLowEdge(iBin+1);
      break;
    }
  }

  return std::make_pair(left, right);
}

void Helper::CD(TFile* file, const std::string& dirName) {
  if(file == nullptr) throw std::runtime_error("Helper::CD() - file is nullptr");

  if(file->GetDirectory(dirName.c_str()) == nullptr) file->mkdir(dirName.c_str());
  file->cd(dirName.c_str());
}

TGraph* Helper::EvaluateMovingAverage(const TGraph* graphIn, int aveLength, bool excludeOwnPoint) {
  TGraph* graphOut = new TGraph();
  EvaluateMovingAverage(graphIn, graphOut, aveLength);

  return graphOut;
}

void Helper::EvaluateMovingAverage(const TGraph* graphIn, TGraph* graphOut, int radius, bool isExcludeOwnPoint) {
  if(graphIn == nullptr || graphOut == nullptr) throw std::runtime_error("Helper::EvaluateMovingAverage(): graphIn == nullptr || graphOut == nullptr");

  const int nPoints = graphIn->GetN();
  for(int iPoint=0; iPoint<nPoints; ++iPoint) {
    int nLocalPoints{0};
    double value{0.};
    const int localRadius = std::min(radius, std::min(iPoint, nPoints-1-iPoint));
    for(int iLocalPoint=iPoint-localRadius; iLocalPoint<=iPoint+localRadius; ++iLocalPoint) {
      if(isExcludeOwnPoint && iLocalPoint == iPoint && localRadius!=0) continue;
      value += graphIn->GetPointY(iLocalPoint);
      ++nLocalPoints;
    }
    graphOut->SetPoint(iPoint, graphIn->GetPointX(iPoint), value/nLocalPoints);
  } // nPoints
}

void Helper::DivideGraph(TGraph* num, const TGraph* den) {
  if(num == nullptr || den == nullptr) throw std::runtime_error("Helper::DivideGraph(): num == nullptr || den == nullptr");
  const int nPoints = num->GetN();
  if(den->GetN() != nPoints) throw std::runtime_error("Helper::DivideGraph(): den->GetN() != nPoints");
  for(int iPoint=0; iPoint<nPoints; ++iPoint) {
    const double xNum = num->GetPointX(iPoint);
    const double xDen = den->GetPointX(iPoint);
    const double yNum = num->GetPointY(iPoint);
    const double yDen = den->GetPointY(iPoint);
    if(std::fabs(xNum - xDen)>1e-4) throw std::runtime_error("Helper::DivideGraph(): num and den points X coordinates do not match");
    num->SetPointY(iPoint, yNum/yDen);
  }
}

TF1* Helper::FitLifetimeHisto(TH1* histo, const std::string& option) {
  const double lo = histo->GetBinLowEdge(1) + 1e-3;
  const double hi = histo->GetBinLowEdge(histo->GetNbinsX()+1) - 1e-3;
  auto parEst = Helper::EstimateExpoParameters(histo, lo, hi);
  TF1* fitFunc = new TF1("fitFunc", "[0]*TMath::Exp(-x/[1])", lo, hi);
  fitFunc->SetParameters(parEst.first, parEst.second);
  histo->Fit(fitFunc, ("0"+option).c_str(), "", lo, hi);
  fitFunc->SetLineColor(histo->GetLineColor());

  return fitFunc;
}

void Helper::DivideFunctionByHisto(TH1* histo, TF1* func, const std::string& option) {
  if(option.empty()) {
    histo->Divide(func);
  } else if(option == "I") {
    for(int iBin=1, nBins=histo->GetNbinsX(); iBin<=nBins; iBin++) {
      const double histoValue = histo->GetBinContent(iBin);
      const double histoError = histo->GetBinError(iBin);
      const double lo = histo->GetBinLowEdge(iBin);
      const double hi = histo->GetBinLowEdge(iBin+1);
      const double funcAverage = func->Integral(lo, hi) / (hi-lo);
      histo->SetBinContent(iBin, histoValue/funcAverage);
      histo->SetBinError(iBin, histoError/funcAverage);
    }
  } else {
    throw std::runtime_error("Helper::DivideFunctionByHisto() - 'option' must be either empty string or I");
  }
}