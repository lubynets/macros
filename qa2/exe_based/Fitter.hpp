//
// Created by oleksii on 04.02.25.
//

#ifndef QA2_FITTER_HPP
#define QA2_FITTER_HPP

#include <TH1.h>

namespace Fitter {

void Fit(TH1* h);

void FitPeak(TH1* h);

//======================================================================================================================

void Fit(TH1* h) {
  FitPeak(h);
}

void FitPeak(TH1* h) {
  TF1* func = new TF1("func", "gaus", )
}

}

#endif //QA2_FITTER_HPP
