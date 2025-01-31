#include "Helper.hpp"

#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TGraphMultiErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TPaveText.h>
#include <TROOT.h>

#include <vector>

struct HistoQuantities {
    float underflow_{-999.f};
    float overflow_{-999.f};
    float mean_{-999.f};
    float mean_err_{-999.f};
    float stddev_{-999.f};
    float stddev_err_{-999.f};
};

std::vector<std::pair<std::string, std::string>> FindCuts(const TFile* fileIn, std::string name_start);
bool stofCompare(std::pair<std::string, std::string> a, std::pair<std::string, std::string> b);
void SetLineDrawParameters(std::vector<TF1*> fs, int lineWidth=1, int lineStyle=7, Color_t lineColor=kBlack);
void CustomizeGraphYRange(TGraphMultiErrors* graph);
HistoQuantities EvaluateHistoQuantities(const TH1* h);
TPaveText ConvertHistoQuantitiesToText(const HistoQuantities& q, float x1, float y1, float x2, float y2);

void mc_qa2diff(const std::string& fileName, int prompt_or_nonprompt) {
  TString currentMacroPath = __FILE__;
  TString directory = currentMacroPath(0, currentMacroPath.Last('/'));
  gROOT->Macro( directory + "/mc_qa2.style.cc" );

  if(prompt_or_nonprompt !=1 && prompt_or_nonprompt != 2) {
    throw std::runtime_error("prompt_or_nonprompt must be 1 or 2");
  }

  const std::string promptness = prompt_or_nonprompt == 1 ? "prompt" : "nonprompt";

  TFile* fileIn = TFile::Open(fileName.c_str());
  if(fileIn == nullptr) {
    throw std::runtime_error("fileIn == nullptr");
  }

  std::string fileOutName = "mc_qa2diff";

  struct Variable {
    std::string name_;
    std::string cut_name_;
    std::string cut_title_;
    std::string cut_unit_;
    bool log_res_;
    bool log_pull_;
  };

  std::vector<Variable> vars {
//  name    cutname  title          unit    logres logpull
    {"P",   "psim",  "p^{mc}",     "GeV/c", false, false},
    {"Pt",  "pTsim", "p_{T}^{mc}", "GeV/c", false, false},
    {"Xsv", "psim",  "p^{mc}",     "GeV/c", false, false},
    {"Ysv", "psim",  "p^{mc}",     "GeV/c", false, false},
    {"Zsv", "psim",  "p^{mc}",     "GeV/c", false, false},
    {"L",   "lsim",  "L^{mc}",     "cm",    false, false},
    {"T",   "tsim",  "T^{mc}",     "ps",    false, false},
  };

  struct ResPull{
    std::string name_;
    std::string prefix_;
    bool is_draw_oneline_on_stddev_;
  };

  std::vector<ResPull> resPulls {
    {"Residual", "res",  false},
    {"Pull",     "pull", true}
  };

  for(auto& var : vars) {
    bool is_first_canvas{true};
    std::string printing_bracket = "(";

    auto cutsVar = FindCuts(fileIn, "Candidates_Simulated_"  + promptness + "_" + var.cut_name_);

    std::vector<TGraphMultiErrors*> grMu;
    std::vector<TGraphErrors*> grSigma;
    grMu.resize(resPulls.size());
    grSigma.resize(resPulls.size());
    for(int iRP=0; iRP<resPulls.size(); iRP++) {
      grMu.at(iRP) = new TGraphMultiErrors(cutsVar.size(), 2);
      grSigma.at(iRP) = new TGraphErrors(cutsVar.size());
      for(auto& g : (std::vector<TGraph*>){grMu.at(iRP), grSigma.at(iRP)}) {
        g->SetTitle("");
        g->GetXaxis()->SetTitle((var.cut_title_ + ", " + var.cut_unit_).c_str());
        g->SetMarkerStyle(kFullSquare);
        g->SetMarkerSize(2);
        g->SetMarkerColor(kBlue);
        g->SetLineWidth(3);
        g->SetLineColor(kBlue);
        g->SetFillColorAlpha(kBlue, 0.3);
      }
    }

    TPaveText generalText(0.74, 0.82, 0.87, 0.90, "brNDC");
    generalText.SetFillColor(0);
    generalText.SetTextSize(0.03);
    generalText.SetTextFont(62);
    generalText.AddText(promptness.c_str());

    int iPoint{0};
    for(auto& cV : cutsVar) {
      TPaveText cutText = generalText;
      cutText.AddText((cV.first + " < " + var.cut_title_ + " < " + cV.second + " " + var.cut_unit_).c_str());
      const float grX = (atof(cV.first.c_str()) + atof(cV.second.c_str())) / 2;
      const float grEX = -(atof(cutsVar.at(0).second.c_str()) - atof(cutsVar.at(0).first.c_str())) / 15;

      for(int iRP=0; iRP<resPulls.size(); iRP++) {
        TH1D* hIn = fileIn->Get<TH1D>(("Candidates_Simulated_"  + promptness + "_"  + var.cut_name_ + "_"  + cV.first + "_"  + cV.second + "/" +  resPulls.at(iRP).prefix_ + "_"  + var.name_ + "_"  + var.cut_name_ + "_"  + cV.first + "_"  + cV.second).c_str());
        if(hIn == nullptr) {
          throw std::runtime_error("hIn == nullptr for " + var.name_ + "; " + var.cut_name_ + "; " + cV.first + "; "  + cV.second + "; " + resPulls.at(iRP).prefix_);
        }

        if(is_first_canvas && iPoint==0) {
          grMu.at(iRP)->GetYaxis()->SetTitle((resPulls.at(iRP).name_ + " mean of " + hIn->GetXaxis()->GetTitle()).c_str());
          grSigma.at(iRP)->GetYaxis()->SetTitle((resPulls.at(iRP).name_ + " width of " + hIn->GetXaxis()->GetTitle()).c_str());
        }

        hIn->UseCurrentStyle();

        TCanvas cc(("cc" + resPulls.at(iRP).prefix_).c_str(), ("cc" + resPulls.at(iRP).prefix_).c_str(), 1200, 800);
        cc.SetLogy(var.log_res_);
        hIn->Draw();
        cutText.Draw("same");
        HistoQuantities quant = EvaluateHistoQuantities(hIn);
        TPaveText quant_text = ConvertHistoQuantitiesToText(quant, 0.74, 0.6, 0.87, 0.7);
        quant_text.Draw("same");

        grMu.at(iRP)->SetPoint(iPoint, grX, quant.mean_);
        grMu.at(iRP)->SetPointEX(iPoint, grEX, grEX);
        grMu.at(iRP)->SetPointEY(iPoint, 0, quant.mean_err_, quant.mean_err_);
        grMu.at(iRP)->SetPointEY(iPoint, 1, quant.stddev_, quant.stddev_);

        grSigma.at(iRP)->SetPoint(iPoint, grX, quant.stddev_);
        grSigma.at(iRP)->SetPointError(iPoint, 0, quant.stddev_err_);

        if(is_first_canvas) printing_bracket = "(";
        else                printing_bracket = "";
        cc.Print((fileOutName + "_" + var.name_ + "_" +  resPulls.at(iRP).prefix_ + ".pdf" + printing_bracket).c_str(), "pdf");
      } // resPulls
      is_first_canvas = false;
      ++iPoint;
    } // cutsVar

    const float xlo = grMu.at(0)->GetPointX(0);
    const float xhi = grMu.at(0)->GetPointX(grMu.at(0)->GetN()-1);
    TF1 zeroLine("zeroLine", "[0]", xlo, xhi);
    TF1 oneLine("oneLine", "[0]", xlo, xhi);
    zeroLine.SetParameter(0, 0);
    oneLine.SetParameter(0, 1);
    SetLineDrawParameters({&zeroLine, &oneLine});

    TCanvas emptycanvas("emptycanvas", "", 1200, 800);
    for(int iRP=0; iRP<resPulls.size(); iRP++) {
      TLegend legStat(0.24, 0.83, 0.37, 0.90);
      TLegend legBoth(0.24, 0.71, 0.37, 0.90);
      legStat.AddEntry(grMu.at(0), "stat. error", "E");
      legBoth.AddEntry(grMu.at(0), "stat. error", "E");
      auto entry = legBoth.AddEntry("", "distrib. width", "F");
      entry->SetFillColorAlpha(grMu.at(0)->GetFillColor(), 0.3);
      entry->SetLineColor(kWhite);
      entry->SetFillStyle(1000);

      TCanvas ccMuStat(("cc" +  resPulls.at(iRP).prefix_ + "Mu").c_str(), ("cc" +  resPulls.at(iRP).prefix_ + "Mu").c_str(), 1200, 800);
      grMu.at(iRP)->Draw("AP; 2");
      zeroLine.Draw("L same");
      legBoth.Draw("same");
      ccMuStat.Print((fileOutName + "_" + var.name_ + "_" +  resPulls.at(iRP).prefix_ + ".pdf").c_str(), "pdf");

      TCanvas ccMuWidth(("cc" +  resPulls.at(iRP).prefix_ + "Mu").c_str(), ("cc" +  resPulls.at(iRP).prefix_ + "Mu").c_str(), 1200, 800);
      grMu.at(iRP)->Draw("AP; X");
      CustomizeGraphYRange(grMu.at(iRP));
      zeroLine.Draw("L same");
      legStat.Draw("same");
      ccMuWidth.Print((fileOutName + "_" + var.name_ + "_" +  resPulls.at(iRP).prefix_ + ".pdf").c_str(), "pdf");

      TCanvas ccSigma(("cc" +  resPulls.at(iRP).prefix_ + "Sigma").c_str(), ("cc" +  resPulls.at(iRP).prefix_ + "Sigma").c_str(), 1200, 800);
      grSigma.at(iRP)->Draw("AP");
      if(resPulls.at(iRP).is_draw_oneline_on_stddev_) oneLine.Draw("L same");
      legStat.Draw("same");
      ccSigma.Print((fileOutName + "_" + var.name_ + "_" +  resPulls.at(iRP).prefix_ + ".pdf").c_str(), "pdf");

      emptycanvas.Print((fileOutName + "_" + var.name_ + "_" +  resPulls.at(iRP).prefix_ + ".pdf]").c_str(), "pdf");
    } // resPulls
  } // vars

  fileIn->Close();
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./mc_qa2diff fileName (prompt_or_nonprompt)" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileName = argv[1];
  const int prompt_or_nonprompt = argc>2 ? atoi(argv[2]) : 1;

  mc_qa2diff(fileName, prompt_or_nonprompt);

  return 0;
}

HistoQuantities EvaluateHistoQuantities(const TH1* h) {
  HistoQuantities result;
  const float integral = h->GetEntries();
  result.underflow_ = h->GetBinContent(0) / integral;
  result.overflow_ = h->GetBinContent(h->GetNbinsX()+1) / integral;
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

  text.AddText(("underflow = " + Helper::to_string_with_precision(q.underflow_*100, 2) + "%").c_str());
  text.AddText(("overflow = " + Helper::to_string_with_precision(q.overflow_*100, 2) + "%").c_str());
  text.AddText(("#mu = " + Helper::to_string_with_precision(q.mean_, 3) + " #pm " + Helper::to_string_with_precision(q.mean_err_, 3) + " (stat.)").c_str());
  text.AddText(("#sigma = " + Helper::to_string_with_precision(q.stddev_, 3) + " #pm " + Helper::to_string_with_precision(q.stddev_err_, 3) + " (stat.)").c_str());

  return text;
}

std::vector<std::pair<std::string, std::string>> FindCuts(const TFile* fileIn, std::string name_start) {
  if(name_start.back() != '_') name_start.push_back('_');
  std::vector<std::pair<std::string, std::string>> result;

  auto lok = fileIn->GetListOfKeys();
  const int nDirs = lok->GetEntries();
  for(int iDir=0; iDir<nDirs; iDir++) {
    const std::string dirName = lok->At(iDir)->GetName();
    if(dirName.substr(0, name_start.size()) != name_start) continue;
    std::pair<std::string, std::string> cutPair;
    bool isFirstCutRead{false};
    for(int iChar=name_start.size(); iChar<dirName.size(); iChar++) {
      char letter = dirName.at(iChar);
      if(letter != '_') {
        if(!isFirstCutRead) cutPair.first.push_back(letter);
        else                cutPair.second.push_back(letter);
      } else {
        isFirstCutRead = true;
      }
    }
    result.emplace_back(cutPair);
  }

  std::sort(result.begin(), result.end(), stofCompare);

  if(result.size() == 0) {
    throw std::runtime_error("FindCuts(): " + name_start + " cuts are not present");
  }

  return result;
}

bool stofCompare(std::pair<std::string, std::string> a, std::pair<std::string, std::string> b) {
  return atof(a.first.c_str()) < atof(b.first.c_str());
}

void SetLineDrawParameters(std::vector<TF1*> fs, int lineWidth, int lineStyle, Color_t lineColor) {
  for(auto& f : fs) {
    f->SetLineWidth(lineWidth);
    f->SetLineStyle(lineStyle);
    f->SetLineColor(lineColor);
  }
}

void CustomizeGraphYRange(TGraphMultiErrors* graph) {
  const int nPoints = graph->GetN();
  float min = 1e9;
  float max = -1e9;

  for(int iPoint=0; iPoint<nPoints; iPoint++) {
    const float up = graph->GetPointY(iPoint) + graph->GetErrorY(iPoint, 0);//only 1-st error is needed
    const float lo = graph->GetPointY(iPoint) - graph->GetErrorY(iPoint, 0);//only 1-st error is needed

    min = std::min(min, lo);
    max = std::max(max, up);
  }

  const float diff = max-min;
  max += diff/10;
  min -= diff/10;

  graph->GetYaxis()->SetRangeUser(min, max);
}
