//
// Created by oleksii on 13.10.2026
//

#ifndef QA2_FORWARD_FOLDED_EXPO_HPP
#define QA2_FORWARD_FOLDED_EXPO_HPP

#include <TH2.h>

class ForwardFoldedExpo {
public:
  void SetResponseMatrix(TH2* rm) { response_matrix_ = rm; }
  void Init();

  double operator()(double* x, double* par);

private:
  double IntegrateExpo(double from, double to);
  double IntegrateExpo(int iBin);

  int n_bins_{};
  double expo_par_A_{};
  double expo_par_tau_{};
  TH2* response_matrix_{};
};

#endif //QA2_FORWARD_FOLDED_EXPO_HPP
