//
// Created by oleksii on 05.02.25.
//

#ifndef QA2_SHAPES_HPP
#define QA2_SHAPES_HPP

#include <TMath.h>

namespace PolN {
enum Pars : int {
  kShift = 0,
  kConst,
  kSlope,
  kSecond,
  kThird,
  nPars
};

inline double Shape(const double* x, const double* par) {
  double shift = par[kShift]; // designed to be fixed at expected peak position
  double a0 = par[kConst];
  double a1 = par[kSlope];
  double a2 = par[kSecond];
  double a3 = par[kThird];

  double u = x[0] - shift;

  return a0 + a1 * u + a2 * u*u + a3 * u*u*u;
}
}

namespace Gaus {
enum Pars : int {
  kFactor = 0,
  kShift,
  kMu,
  kSigma,
  nPars
};

inline double Shape(const double* x, const double* par) {
  double factor = par[kFactor];
  double shift = par[kShift]; // designed to be fixed at expected peak position
  double mu = par[kMu];
  double sigma = par[kSigma];

  double u = (x[0] - shift - mu) / sigma;

  return factor * TMath::Exp(-u * u / 2);
}
}

namespace DoubleGaus {
enum Pars : int {
  kFactor1 = 0,
  kFactor2,
  kShift,
  kMu,
  kSigma1,
  kSigma2,
  nPars
};

inline double Shape(const double* x, const double* par) {
  double factor1 = par[kFactor1];
  double factor2 = par[kFactor2];
  double shift = par[kShift]; // designed to be fixed at expected peak position
  double mu = par[kMu];
  double sigma1 = par[kSigma1];
  double sigma2 = par[kSigma2];

  double u1 = (x[0] - shift - mu) / sigma1;
  double u2 = (x[0] - shift - mu) / sigma2;

  return factor1 * TMath::Exp(-u1 * u1 / 2) + factor2 * TMath::Exp(-u2 * u2 / 2);
}
}

namespace DoubleSidedCrystalBall {
enum Pars : int {
  kFactor = 0,
  kShift,
  kMu,
  kSigma,
  kA1,
  kN1,
  kA2,
  kN2,
  nPars
};

inline double Shape(const double* x, const double* par) {
  double factor = par[kFactor];
  double shift = par[kShift]; // designed to be fixed at expected peak position
  double mu = par[kMu];
  double sigma = par[kSigma];
  double a1 = par[kA1];
  double n1 = TMath::Power(10, par[kN1]);
  double a2 = par[kA2];
  double n2 = TMath::Power(10, par[kN2]);

  double u = (x[0] - shift - mu) / sigma;

  if (u < -a1)
    return factor * TMath::Exp(-a1 * a1 / 2) * TMath::Power(1 - a1 * (u + a1) / n1, -n1);
  else if (u >= -a1 && u < a2)
    return factor * TMath::Exp(-u * u / 2);
  else if (u >= a2)
    return factor * TMath::Exp(-a2 * a2 / 2) * TMath::Power(1 + a2 * (u - a2) / n2, -n2);
  else
    return -1.;
}
}

#endif //QA2_SHAPES_HPP
