//
// Created by oleksii on 21.07.25.
//

#include "HelperPlot.hpp"

#include "HelperGeneral.hpp"

#include <TF1.h>

HelperPlot::HistoQuantities HelperPlot::EvaluateHistoQuantities(const TH1* h) {
  HelperPlot::HistoQuantities result;
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

TPaveText HelperPlot::ConvertHistoQuantitiesToText(const HelperPlot::HistoQuantities& q, float x1, float y1, float x2, float y2) {
  TPaveText text(x1, y1, x2, y2, "brNDC");
  text.SetFillColor(0);
  text.SetTextSize(0.03);
  text.SetTextFont(62);

  text.AddText(("nentries = " + HelperGeneral::to_string_with_precision(q.nentries_, 0)).c_str());
  text.AddText(("underflow = " + HelperGeneral::to_string_with_precision(q.underflow_ * 100, 2) + "%").c_str());
  text.AddText(("overflow = " + HelperGeneral::to_string_with_precision(q.overflow_ * 100, 2) + "%").c_str());
  text.AddText(
          ("#mu = " + HelperGeneral::to_string_with_precision(q.mean_, 3) + " #pm " + HelperGeneral::to_string_with_precision(q.mean_err_, 3) +
           " (stat.)").c_str());
  text.AddText(("#sigma = " + HelperGeneral::to_string_with_precision(q.stddev_, 3) + " #pm " +
                HelperGeneral::to_string_with_precision(q.stddev_err_, 3) + " (stat.)").c_str());

  return text;
}

void HelperPlot::CustomizeGraphYRange(TGraphMultiErrors* graph, int ne, TF1* f) {
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

void HelperPlot::SetLineDrawParameters(std::vector<TF1*> fs, int lineWidth, int lineStyle, Color_t lineColor) {
  for (auto& f: fs) {
    f->SetLineWidth(lineWidth);
    f->SetLineStyle(lineStyle);
    f->SetLineColor(lineColor);
  }
}

TF1* HelperPlot::HorizontalLine4Graph(float level, TGraph* graph) {
  float xlo = graph->GetPointX(0);
  float xhi = graph->GetPointX(graph->GetN() - 1);
  const float diff = xhi - xlo;
  xlo = -diff / 10;
  xhi += diff / 10;
  TF1* horizLine = new TF1("horizLine", "[0]", xlo, xhi);
  horizLine->SetParameter(0, level);

  return horizLine;
}

TPaveText* HelperPlot::AddOneLineText(const std::string& text, const std::array<float, 4>& xy, const std::string& option, float size) {
  if(text.empty()) return nullptr;
  TPaveText* textPtr = new TPaveText(xy.at(0), xy.at(1), xy.at(2), xy.at(3), option.c_str());
  textPtr->SetFillColor(0);
  textPtr->SetTextSize(size);
  textPtr->SetTextFont(62);
  textPtr->AddText(text.c_str());
  textPtr->Draw("same");
  return textPtr;
}

std::vector<TPaveText*> HelperPlot::AddMultiLineText(const std::vector<std::string>& texts, const std::array<float, 4>& xy, const std::string& option, float size) {
  const float step = xy.at(3) - xy.at(1);
  std::vector<TPaveText*> ptr;
  for(int iLine=0, nLines=texts.size(); iLine<nLines; iLine++) {
    ptr.emplace_back(HelperPlot::AddOneLineText(texts.at(iLine), {xy.at(0), xy.at(1)-iLine*step, xy.at(2), xy.at(3)-iLine*step}, option, size));
  }
  return ptr;
}

std::pair<double, double> HelperPlot::GetMinMaxBinWithError(const TH1* h) {
  double min = 1e9;
  double max = -1e9;
  for(int iBin=1; iBin<=h->GetNbinsX(); iBin++) {
    min = std::min(min, h->GetBinContent(iBin) - h->GetBinError(iBin));
    max = std::max(max, h->GetBinContent(iBin) + h->GetBinError(iBin));
  }
  return std::make_pair(min, max);
}

void HelperPlot::SlightlyShiftXAxis(TGraph* gr, float value) {
  if(value == -1) {
    value = (gr->GetPointX(gr->GetN()-1) - gr->GetPointX(0)) / 50;
  }
  const int nPoints = gr->GetN();
  for(int iPoint=0; iPoint<nPoints; iPoint++) {
    gr->SetPointX(iPoint, gr->GetPointX(iPoint) + value);
  }
}

void HelperPlot::ScaleCanvasVertically(TCanvas* cTo, const TCanvas* cFrom, double scaleFactor) {
  cTo->SetCanvasSize(cFrom->GetWw(), cFrom->GetWh()/scaleFactor);
  cTo->SetBottomMargin(cFrom->GetBottomMargin()*scaleFactor);
  cTo->SetTopMargin(cFrom->GetTopMargin()*scaleFactor);
}

void HelperPlot::CustomizeHistogramsYRange(const std::vector<TH1*>& histos, bool isLog, double lo, double hi, double part) {
  double max = -1e9;
  double min = 1e9;
  for(auto& histo : histos) {
    auto [hmin, hmax] = HelperPlot::GetMinMaxBinWithError(histo);
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

void HelperPlot::CloseCanvasPrinting(const std::vector<std::string>& names) {
  TCanvas emptyCanvas("emptyCanvas", "", 1000, 1000);
  for(auto& name : names) {
    emptyCanvas.Print((name + ".pdf]").c_str(), "pdf");
  }
}