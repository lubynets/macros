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
constexpr double massLambdaCDetectorWidth{0.01}; // by order of magnitude. Needed for initializing the fit.

constexpr double UndefValueDouble{-999.};
constexpr float UndefValueFloat{-999.f};
constexpr int UndefValueInt{-999};

struct HistoQuantities {
  float nentries_{-999.f};
  float underflow_{-999.f};
  float overflow_{-999.f};
  float mean_{-999.f};
  float mean_err_{-999.f};
  float stddev_{-999.f};
  float stddev_err_{-999.f};
};

std::string getSubstringBeforeLastSlash(const std::string& input);

std::vector<std::pair<std::string, std::string>> FindCuts(TFile* fileIn, std::string name_start, bool printCuts=false);

bool stofCompare(std::pair<std::string, std::string> a, std::pair<std::string, std::string> b);

HistoQuantities EvaluateHistoQuantities(const TH1* h);

TPaveText ConvertHistoQuantitiesToText(const HistoQuantities& q, float x1, float y1, float x2, float y2);

void CustomizeGraphYRange(TGraphMultiErrors* graph, int ne = 1, TF1* f = nullptr);

void CustomizeHistogramsYRange(const std::vector<TH1*>& histos, double lo=-1e9, double hi=1e9, double part = 0.9);

std::pair<double, double> GetMinMaxBinWithError(const TH1* h);

void SetLineDrawParameters(std::vector<TF1*> fs, int lineWidth = 1, int lineStyle = 7, Color_t lineColor = kBlack);

TF1* HorizontalLine4Graph(float level, TGraph* graph);

bool string_to_bool(const std::string& str);

TFile* OpenFileWithNullptrCheck(const std::string& fileName, const std::string& option = "read");

template<typename T>
T* GetObjectWithNullptrCheck(TFile* fileIn, const std::string& objectName) {
  T* ptr = fileIn->Get<T>(objectName.c_str());
  if(ptr == nullptr) {
    throw std::runtime_error("GetObjectWithNullptrCheck() - object " + objectName + " in file " + fileIn->GetName() + " is missing");
  }
  return ptr;
}

template<typename T>
inline std::string to_string_with_precision(const T a_value, const int n=2) {
  std::ostringstream out;
  out.precision(n);
  out << std::fixed << a_value;
  return out.str();
}

template<typename T>
inline std::string to_string_with_significant_figures(const T a_value, const int n=2) {
  if(a_value == 0) return "0";

  const double dMag = std::log10(std::abs(a_value)); // scale of the a_value (e.g 1.* for 1.2345, 2.* for 12.345 etc)
  const int iMag = static_cast<int>(dMag-n+1 > 0 ? dMag-n+1 : dMag-n);
  const T shifted_value = a_value/std::pow(10, iMag); // shift decimal point to have all required digits to l.h.s. from it
  const T rounded_value = static_cast<T>(std::round(shifted_value)); // get rid of r.h.s. from decimal point
  const T reshifted_value = rounded_value*std::pow(10, iMag); // return decimal point to its original place
  const int precision = iMag < 0 ? -iMag : 0; // determine how many digits after decimal point one needs
  return to_string_with_precision(reshifted_value, precision);
}

TPaveText* AddOneLineText(const std::string& text, const std::array<float, 4>& xy, const std::string& option="brNDC", float size=0.03);

void SlightlyShiftXAxis(TGraph* gr, float value = -1);

template<typename T>
inline void RemoveEdgeLabelFromAxis(T* obj, const std::string& edge, const std::string& axisletter) {
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

std::pair<float, float> EstimateExpoParameters(TH1* h, float lo, float hi);

void PrintInfoOnTF1(const TF1* f);

void CloseCanvasPrinting(const std::vector<std::string>& names);

void CheckHistogramsForXaxisIdentity(const TH1* h1, const TH1* h2);

} // namespace Helper
#endif //QA2_HELPER_HPP
