//
// Created by oleksii on 11.02.25.
//
#include "Helper.hpp"

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

std::vector<std::pair<std::string, std::string>> Helper::FindCuts(TFile* fileIn, std::string name_start) {
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
    for (int iChar = name_start.size(); iChar < fullDirName.size(); iChar++) {
      char letter = fullDirName.at(iChar);
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

  if (result.empty()) {
    throw std::runtime_error("FindCuts(): " + name_start + " cuts are not present");
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

void Helper::AddOneLineText(const std::string& text, float x1, float y1, float x2, float y2) {
  if(text.empty()) return;
  TPaveText* textPtr = new TPaveText(x1, y1, x2, y2, "brNDC");
  textPtr->SetFillColor(0);
  textPtr->SetTextSize(0.03);
  textPtr->SetTextFont(62);
  textPtr->AddText(text.c_str());
  textPtr->Draw("same");
}

bool Helper::string_to_bool(const std::string& str) {
  if(str == "true") return true;
  else if(str == "false") return false;
  else throw std::runtime_error("string_to_bool(): argument must be either true or false");
}

void Helper::CustomizeHistogramsYRange(const std::vector<TH1*>& histos, double absoluteMaximum) {
  double max = -1e9;
  for(auto& histo : histos) {
    max = std::max(max, histo->GetMaximum());
  }
  max *= 1.1;
  max = std::min(max, absoluteMaximum);
  for(auto& histo : histos) {
    histo->GetYaxis()->SetRangeUser(histo->GetMinimum(), max);
  }
}