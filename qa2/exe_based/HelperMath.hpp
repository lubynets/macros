//
// Created by oleksii on 21.07.25.
//

#ifndef QA2_HELPERMATH_HPP
#define QA2_HELPERMATH_HPP

#include <TGraph.h>

namespace HelperMath {

std::pair<float, float> EstimateExpoParameters(TH1* h, float lo, float hi);

std::pair<double, double> DetermineWorkingRangesTH1(const TH1* histo, double leftMargin=0.0015, double rightMargin=0.0015);

TGraph* EvaluateMovingAverage(const TGraph* graphIn, int radius, bool isExcludeOwnPoint=false);
void EvaluateMovingAverage(const TGraph* graphIn, TGraph* graphOut, int radius, bool isExcludeOwnPoint=false);

void DivideGraph(TGraph* num, const TGraph* den);

template <typename T, size_t N>
struct nested_vector {
    using type = std::vector<typename nested_vector<T, N - 1>::type>;
};
// Base case
template <typename T>
struct nested_vector<T, 0> {
    using type = T;
};

template <typename T, std::size_t N>
using tensor = typename nested_vector<T, N>::type;

template <typename T, std::size_t N>
tensor<T, N> make_tensor(std::array<std::size_t, N> const& dims, T const& init = {}) {
  if constexpr (N == 0) {
    return init;
  } else {
  tensor<T, N> v;
  v.resize(dims[0]);
  auto subdims = [&] {
      std::array<std::size_t, N - 1> tail{};
      std::copy(dims.begin() + 1, dims.end(), tail.begin());
      return tail;
  }();
  for (auto& sub : v) {
    sub = make_tensor<T, N - 1>(subdims, init);
  }
  return v;
  }
}

template<typename T>
using tensor2 = tensor<T, 2>;

template<typename T>
using tensor3 = tensor<T, 3>;

TF1* FitLifetimeHisto(TH1* histo, const std::string& option="");

void DivideFunctionByHisto(TH1* histo, TF1* func, const std::string& option="");

std::pair<TH1*, TH1*> EvaluateEfficiencyHisto(TH1* hNum, TH1* hDen);

TH1* MergeHistograms(const std::vector<TH1*>& histos);

};


#endif //QA2_HELPERMATH_HPP
