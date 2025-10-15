//
// Created by oleksii on 21.07.25.
//

#ifndef QA2_HELPERPLOT_HPP
#define QA2_HELPERPLOT_HPP

#include <TCanvas.h>
#include <TGraphMultiErrors.h>
#include <TH1.h>
#include <TPaveText.h>

#include <stdexcept>

namespace HelperPlot {

struct HistoQuantities {
    float nentries_{-999.f};
    float underflow_{-999.f};
    float overflow_{-999.f};
    float mean_{-999.f};
    float mean_err_{-999.f};
    float stddev_{-999.f};
    float stddev_err_{-999.f};
};

HistoQuantities EvaluateHistoQuantities(const TH1* h);

TPaveText ConvertHistoQuantitiesToText(const HistoQuantities& q, float x1, float y1, float x2, float y2);

void CustomizeGraphYRange(TGraphMultiErrors* graph, int ne = 1, TF1* f = nullptr);

void CustomizeHistogramsYRange(const std::vector<TH1*>& histos, bool isLog=false, double lo=-1e9, double hi=1e9, double part = 0.9);

std::pair<double, double> GetMinMaxBinWithError(const TH1* h);

void SetLineDrawParameters(std::vector<TF1*> fs, int lineWidth = 1, int lineStyle = 7, Color_t lineColor = kBlack);

TF1* HorizontalLine4Graph(double level, TGraph* graph);

TPaveText* AddOneLineText(const std::string& text, const std::array<float, 4>& xy, const std::string& option="brNDC", float size=0.03);

std::vector<TPaveText*> AddMultiLineText(const std::vector<std::string>& texts, const std::array<float, 4>& xy, const std::string& option="brNDC", float size=0.03);

void SlightlyShiftXAxis(TGraph* gr, float value = -1);

template<typename T>
inline void RemoveEdgeLabelFromAxis(T* obj, const std::string& edge, const std::string& axisletter) {
  if(axisletter != "x" && axisletter != "y") {
    throw std::runtime_error("HelperPlot::RemoveEdgeLabelFromAxis(): axisletter must be x or y");
  }
  if(edge != "first" && edge != "last") {
    throw std::runtime_error("HelperPlot::RemoveEdgeLabelFromAxis(): edge must be first or last");
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

template<typename T1, typename T2>
inline void ScalePlotVertically(T1* plotTo, const T2* plotFrom, double scaleFactor) {
  plotTo->GetXaxis()->SetTitleSize(plotFrom->GetXaxis()->GetTitleSize()*scaleFactor);
  plotTo->GetXaxis()->SetLabelSize(plotFrom->GetXaxis()->GetLabelSize()*scaleFactor);
  plotTo->GetXaxis()->SetTickLength(plotFrom->GetXaxis()->GetTickLength()*scaleFactor);
  plotTo->GetYaxis()->SetTitleSize(plotFrom->GetYaxis()->GetTitleSize()*scaleFactor);
  plotTo->GetYaxis()->SetLabelSize(plotFrom->GetYaxis()->GetLabelSize()*scaleFactor);
  plotTo->GetYaxis()->SetTitleOffset(plotFrom->GetYaxis()->GetTitleOffset()/scaleFactor);
}

void ScaleCanvasVertically(TCanvas* cTo, const TCanvas* cFrom, double scaleFactor);

void CloseCanvasPrinting(const std::vector<std::string>& names);

inline std::string EvaluatePrintingBracket(size_t vecSize, size_t index) {
  return vecSize == 1 ? "" : index == 0 ? "(" : index == vecSize-1 ? ")" : "";
}

template <typename T>
inline std::string EvaluatePrintingBracket(const std::vector<T>& vec, size_t index) {
  size_t vecSize = vec.size();
  return EvaluatePrintingBracket(vecSize, index);
}

};


#endif //QA2_HELPERPLOT_HPP
