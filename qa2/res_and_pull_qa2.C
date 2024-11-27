template<typename T>
std::string to_string_with_precision(const T a_value, const int n = 6);

struct HistoQuantities {
  float underflow_{-999.f};
  float overflow_{-999.f};
  float mean_{-999.f};
  float mean_err_{-999.f};
  float stddev_{-999.f};
  float stddev_err_{-999.f};
};

HistoQuantities EvaluateHistoQuantities(const TH1* h);
TPaveText ConverHistoQuantitiesToText(const HistoQuantities& q, float x1, float y1, float x2, float y2);

void res_and_pull_qa2(const std::string& fileName, int prompt_or_nonprompt=1) {
  gROOT->Macro( "treeKF_qa2.style.cc" );
  if(prompt_or_nonprompt !=1 && prompt_or_nonprompt != 2) {
    throw std::runtime_error("prompt_or_nonprompt must be 1 or 2");
  }

  const std::string promptness = prompt_or_nonprompt == 1 ? "prompt" : "nonprompt";

  TFile* fileIn = TFile::Open(fileName.c_str());
  if(fileIn == nullptr) {
    throw std::runtime_error("fileIn == nullptr");
  }

  std::string fileOutName = "res_and_pull_qa2";

  struct Variable {
    std::string name_;
  };

  std::vector<Variable> vars {
    {"P" },
    {"Pt"},
    {"X" },
    {"Y" },
    {"Z" },
    {"L" },
    {"T" },
  };

  bool is_first_canvas{true};
  std::string printing_bracket = "(";
  for(auto& var : vars) {
    TH1D* hres = fileIn->Get<TH1D>(("residual/" + promptness + "/hRes_" + var.name_ + "_" + promptness).c_str());
    TH2D* hcorr = fileIn->Get<TH2D>(("corr/" + promptness + "/hCorr_" + var.name_ + "_" + promptness).c_str());
    TH1D* hpull = fileIn->Get<TH1D>(("pull/" + promptness + "/hPull_" + var.name_ + "_" + promptness).c_str());
    if(hres == nullptr || hcorr == nullptr || hpull == nullptr) {
      throw std::runtime_error("hres == nullptr || hcorr == nullptr || hpull == nullptr");
    }

    hres->UseCurrentStyle();
    hcorr->UseCurrentStyle();
    hpull->UseCurrentStyle();

    TPaveText promptnessText(0.14, 0.86, 0.27, 0.90, "brNDC");
    promptnessText.SetFillColor(0);
    promptnessText.SetTextSize(0.03);
    promptnessText.SetTextFont(22);
    promptnessText.AddText(promptness.c_str());

    TCanvas ccRes("ccRes", "ccRes", 1200, 800);
    hres->Draw();
    promptnessText.Draw("same");
    HistoQuantities res_quant = EvaluateHistoQuantities(hres);
    TPaveText res_quant_text = ConverHistoQuantitiesToText(res_quant, 0.14, 0.6, 0.27, 0.7);
    res_quant_text.Draw("same");

    TCanvas ccCorr("ccCorr", "ccCorr", 1200, 800);
    hcorr->Draw();
    promptnessText.Draw("same");

    TCanvas ccPull("ccPull", "ccPull", 1200, 800);
    hpull->Draw();
    promptnessText.Draw("same");
    HistoQuantities pull_quant = EvaluateHistoQuantities(hpull);
    TPaveText pull_quant_text = ConverHistoQuantitiesToText(pull_quant, 0.14, 0.6, 0.27, 0.7);
    pull_quant_text.Draw("same");

    if(is_first_canvas) printing_bracket = "(";
    else                printing_bracket = "";
    ccRes.Print((fileOutName + "_res.pdf" + printing_bracket).c_str(), "pdf");
    ccCorr.Print((fileOutName + "_corr.pdf" + printing_bracket).c_str(), "pdf");
    ccPull.Print((fileOutName + "_pull.pdf" + printing_bracket).c_str(), "pdf");
    is_first_canvas = false;
  }

  printing_bracket = "]";
  TCanvas emptycanvas("emptycanvas", "", 1200, 800);
  emptycanvas.Print((fileOutName + "_res.pdf" + printing_bracket).c_str(), "pdf");
  emptycanvas.Print((fileOutName + "_corr.pdf" + printing_bracket).c_str(), "pdf");
  emptycanvas.Print((fileOutName + "_pull.pdf" + printing_bracket).c_str(), "pdf");
}

template<typename T>
  std::string to_string_with_precision(const T a_value, const int n) {
  std::ostringstream out;
  out.precision(n);
  out << std::fixed << a_value;
  return out.str();
}

HistoQuantities EvaluateHistoQuantities(const TH1* h) {
  HistoQuantities result;
  const float integral = h->GetEntries();
  result.underflow_ = h->GetBinContent(0) / integral;
  result.overflow_ = h->GetBinContent(h->GetNbinsX()) / integral;
  result.mean_ = h->GetMean();
  result.mean_err_ = h->GetMeanError();
  result.stddev_ = h->GetStdDev();
  result.stddev_err_ = h->GetStdDevError();

  return result;
}

TPaveText ConverHistoQuantitiesToText(const HistoQuantities& q, float x1, float y1, float x2, float y2) {
  TPaveText text(x1, y1, x2, y2, "brNDC");
  text.SetFillColor(0);
  text.SetTextSize(0.03);
  text.SetTextFont(22);

  text.AddText(("underflow = " + to_string_with_precision(q.underflow_*100, 2) + "%").c_str());
  text.AddText(("overflow = " + to_string_with_precision(q.overflow_*100, 2) + "%").c_str());
  text.AddText(("#mu = " + to_string_with_precision(q.mean_, 3) + " #pm " + to_string_with_precision(q.mean_err_, 3)).c_str());
  text.AddText(("#sigma = " + to_string_with_precision(q.stddev_, 3) + " #pm " + to_string_with_precision(q.stddev_err_, 3)).c_str());

  return text;
}
