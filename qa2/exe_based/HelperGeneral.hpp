//
// Created by oleksii on 21.07.25.
//

#ifndef QA2_HELPERGENERAL_HPP
#define QA2_HELPERGENERAL_HPP

#include <THnSparse.h>
#include <TFile.h>
#include <TF1.h>

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

namespace HelperGeneral {

constexpr double massLambdaC{2.28646};
constexpr double massLambdaCDetectorWidth{0.01}; // by order of magnitude. Needed for initializing the fit.

constexpr double UndefValueDouble{-999.};
constexpr float UndefValueFloat{-999.f};
constexpr int UndefValueInt{-999};

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

std::vector<std::pair<std::string, std::string>> FindCuts(TFile* fileIn, std::string name_start, bool printCuts=false);

bool string_to_bool(const std::string& str);

TFile* OpenFileWithNullptrCheck(const std::string& fileName, const std::string& option = "read");

template<typename T>
inline T* GetObjectWithNullptrCheck(TFile* fileIn, const std::string& objectName) {
  T* ptr = fileIn->Get<T>(objectName.c_str());
  if(ptr == nullptr) {
    throw std::runtime_error("GetObjectWithNullptrCheck() - object " + objectName + " in file " + fileIn->GetName() + " is missing");
  }
  return ptr;
}

void PrintInfoOnTF1(const TF1* f);

void LoadMacro(const std::string& macroName);

void CD(TFile* file, const std::string& dirName);

void CheckHistogramsForXaxisIdentity(const TH1* h1, const TH1* h2);

std::map<std::string, int> MapAxesIndices(const THnSparse* histo);

void CheckTAxisForRanges(const TAxis& axis, const std::vector<float>& ranges);

void SetTHnSparseAxisRanges(THnSparse* histo, int axisNum, float lo= -999., float hi= -999.);

double InterpolateTH1SuppressWarning(const TH1* h, double value);

void ScaleTHnSparseWithWeight(THnSparse* histoIn, int nDim, const TH1* histoWeight);
};


#endif //QA2_HELPERGENERAL_HPP
