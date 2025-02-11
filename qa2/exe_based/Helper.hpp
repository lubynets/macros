//
// Created by oleksii on 31.01.25.
//

#ifndef QA2_HELPER_HPP
#define QA2_HELPER_HPP

#include <TGraphErrors.h>
#include <TGraphMultiErrors.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TPaveText.h>

#include <string>
#include <sstream>
#include <vector>

namespace Helper {

constexpr double massLambdaC{2.28646};

struct HistoQuantities;

inline std::vector<std::pair<std::string, std::string>> FindCuts(const TFile* fileIn, std::string name_start);

inline bool stofCompare(std::pair<std::string, std::string> a, std::pair<std::string, std::string> b);

inline HistoQuantities EvaluateHistoQuantities(const TH1* h);

inline TPaveText ConvertHistoQuantitiesToText(const HistoQuantities& q, float x1, float y1, float x2, float y2);

inline void CustomizeGraphYRange(TGraphMultiErrors* graph, int ne = 1, TF1* f = nullptr);

inline void SetLineDrawParameters(std::vector<TF1*> fs, int lineWidth = 1, int lineStyle = 7, Color_t lineColor = kBlack);

inline TF1* HorizontalLine4Graph(float level, TGraph* graph);

template<typename T>
inline std::string to_string_with_precision(T a_value, int n = 2);

template<typename T>
inline std::string to_string_with_significant_figures(T a_value, int n=2);

template<typename T>
inline void RemoveEdgeLabelFromAxis(T* obj, const std::string& axisletter);

inline void AddOneLineText(const std::string& text, float x1, float y1, float x2, float y2);

inline void SlightlyShiftXAxis(TGraph* gr, float value = -1);

//======================================================================================================================

void SlightlyShiftXAxis(TGraph* gr, float value) {
  if(value == -1) {
    value = (gr->GetPointX(gr->GetN()-1) - gr->GetPointX(0)) / 50;
  }
  const int nPoints = gr->GetN();
  for(int iPoint=0; iPoint<nPoints; iPoint++) {
    gr->SetPointX(iPoint, gr->GetPointX(iPoint) + value);
  }
}

template<typename T>
std::string to_string_with_precision(const T a_value, const int n) {
  std::ostringstream out;
  out.precision(n);
  out << std::fixed << a_value;
  return out.str();
}

template<typename T>
std::string to_string_with_significant_figures(const T a_value, const int n) {
  if(a_value == 0) return "0";

  const double dMag = std::log10(std::abs(a_value)); // scale of the a_value (e.g 1.* for 1.2345, 2.* for 12.345 etc)
  const int iMag = static_cast<int>(dMag-n+1 > 0 ? dMag-n+1 : dMag-n);
  const T shifted_value = a_value/std::pow(10, iMag); // shift decimal point to have all required digits to l.h.s. from it
  const T rounded_value = static_cast<T>(std::round(shifted_value)); // get rid of r.h.s. from decimal point
  const T reshifted_value = rounded_value*std::pow(10, iMag); // return decimal point to its original place
  const int precision = iMag < 0 ? -iMag : 0; // determine how many digits after decimal point one needs
  return to_string_with_precision(reshifted_value, precision);
}

struct HistoQuantities {
  float nentries_{-999.f};
  float underflow_{-999.f};
  float overflow_{-999.f};
  float mean_{-999.f};
  float mean_err_{-999.f};
  float stddev_{-999.f};
  float stddev_err_{-999.f};
};

HistoQuantities EvaluateHistoQuantities(const TH1* h) {
  HistoQuantities result;
  const float integral = h->GetEntries();
  result.nentries_ = integral;
  result.underflow_ = h->GetBinContent(0) / integral;
  result.overflow_ = h->GetBinContent(h->GetNbinsX() + 1) / integral;
  result.mean_ = h->GetMean();
  result.mean_err_ = h->GetMeanError();
  result.stddev_ = h->GetStdDev();
  result.stddev_err_ = h->GetStdDevError();

  return result;
}

TPaveText ConvertHistoQuantitiesToText(const HistoQuantities& q, float x1, float y1, float x2, float y2) {
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

std::vector<std::pair<std::string, std::string>> FindCuts(const TFile* fileIn, std::string name_start) {
  if (name_start.back() != '_') name_start.push_back('_');
  std::vector<std::pair<std::string, std::string>> result;

  auto lok = fileIn->GetListOfKeys();
  const int nDirs = lok->GetEntries();
  for (int iDir = 0; iDir < nDirs; iDir++) {
    const std::string dirName = lok->At(iDir)->GetName();
    if (dirName.substr(0, name_start.size()) != name_start) continue;
    std::pair<std::string, std::string> cutPair;
    bool isFirstCutRead{false};
    for (int iChar = name_start.size(); iChar < dirName.size(); iChar++) {
      char letter = dirName.at(iChar);
      if (letter != '_') {
        if (!isFirstCutRead) cutPair.first.push_back(letter);
        else cutPair.second.push_back(letter);
      } else {
        isFirstCutRead = true;
      }
    }
    result.emplace_back(cutPair);
  }

  std::sort(result.begin(), result.end(), stofCompare);

  if (result.size() == 0) {
    throw std::runtime_error("FindCuts(): " + name_start + " cuts are not present");
  }

  return result;
}

bool stofCompare(std::pair<std::string, std::string> a, std::pair<std::string, std::string> b) {
  return atof(a.first.c_str()) < atof(b.first.c_str());
}

void CustomizeGraphYRange(TGraphMultiErrors* graph, int ne, TF1* f) {
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

void SetLineDrawParameters(std::vector<TF1*> fs, int lineWidth, int lineStyle, Color_t lineColor) {
  for (auto& f: fs) {
    f->SetLineWidth(lineWidth);
    f->SetLineStyle(lineStyle);
    f->SetLineColor(lineColor);
  }
}

TF1* HorizontalLine4Graph(float level, TGraph* graph) {
  float xlo = graph->GetPointX(0);
  float xhi = graph->GetPointX(graph->GetN() - 1);
  const float diff = xhi - xlo;
  xlo = -diff / 10;
  xhi += diff / 10;
  TF1* horizLine = new TF1("horizLine", "[0]", xlo, xhi);
  horizLine->SetParameter(0, level);

  return horizLine;
}

template<typename T>
void RemoveEdgeLabelFromAxis(T* obj, const std::string& edge, const std::string& axisletter) {
  if(axisletter != "x" && axisletter != "y") {
    throw std::runtime_error("Helper::RemoveEdgeLabelFromAxis(): axisletter must be x or y");
  }
  if(edge != "first" && edge != "last") {
    throw std::runtime_error("Helper::RemoveEdgeLabelFromAxis(): edge must be first or last");
  }
  TAxis* axis = axisletter == "x" ? obj->GetXaxis() : obj->GetYaxis();

  if(axis->IsVariableBinSize()) {
    auto* array = new TArrayD(*axis->GetXbins());
    if(edge == "first") {
      array->SetAt(array->At(0) + 1e-5, 0);
    } else {
      array->SetAt(array->At(array->GetSize()-1) - 1e-5, array->GetSize()-1);
    }
    axis->Set(axis->GetNbins(), array->GetArray());
  } else {
    if(edge == "first") {
      axis->Set(axis->GetNbins(), axis->GetXmin() + 1e-5, axis->GetXmax());
    } else {
      axis->Set(axis->GetNbins(), axis->GetXmin(), axis->GetXmax() - 1e-5);
    }
  }
}

void AddOneLineText(const std::string& text, float x1, float y1, float x2, float y2) {
  TPaveText* textPtr = new TPaveText(x1, y1, x2, y2, "brNDC");
  textPtr->SetFillColor(0);
  textPtr->SetTextSize(0.03);
  textPtr->SetTextFont(62);
  textPtr->AddText(text.c_str());
  textPtr->Draw("same");
}
}
#endif //QA2_HELPER_HPP
