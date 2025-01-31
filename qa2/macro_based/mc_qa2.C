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

void mc_qa2(const std::string& fileName, int prompt_or_nonprompt=1) {
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

  std::string fileOutName = "mc_qa2";

  struct Variable {
    std::string name_;
    bool log_mc_;
    bool log_rec_;
    bool log_res_;
    bool log_corr_;
    bool log_pull_;
  };

  std::vector<Variable> vars {
//  name    logmc  logrec logres logcorr logpull
    {"P",   false, false, false, true, false},
    {"Pt",  false, false, false, true, false},
    {"Xsv", false, false, false, true, false},
    {"Ysv", false, false, false, true, false},
    {"Zsv", false, false, false, true, false},
    {"Xpv", false, false, false, true, false},
    {"Ypv", false, false, false, true, false},
    {"Zpv", false, false, false, true, false},
    {"L",   false, true,  true,  true, true},
    {"T",   false, true,  true,  true, true},
  };

  bool is_first_canvas{true};
  std::string printing_bracket = "(";
  for(auto& var : vars) {
    TH1D* hmc = fileIn->Get<TH1D>(("Simulated_" + promptness + "_total/mc_" + var.name_).c_str());
    TH1D* hrec = fileIn->Get<TH1D>(("Candidates_" + promptness + "_total/rec_" + var.name_).c_str());
    TH1D* hres = fileIn->Get<TH1D>(("Candidates_Simulated_"  + promptness + "_total/res_" + var.name_).c_str());
    TH2D* hcorr = fileIn->Get<TH2D>(("Candidates_Simulated_"  + promptness + "_total/corr_" + var.name_).c_str());
    TH1D* hpull = fileIn->Get<TH1D>(("Candidates_Simulated_"  + promptness + "_total/pull_" + var.name_).c_str());
    if(hmc == nullptr || hrec == nullptr || hres == nullptr || hcorr == nullptr || hpull == nullptr) {
      throw std::runtime_error(("hmc == nullptr || hrec == nullptr || hres == nullptr || hcorr == nullptr || hpull == nullptr for " + var.name_).c_str());
    }

    hmc->UseCurrentStyle();
    hrec->UseCurrentStyle();
    hres->UseCurrentStyle();
    hcorr->UseCurrentStyle();
    hpull->UseCurrentStyle();

    TPaveText promptnessText(0.74, 0.86, 0.87, 0.90, "brNDC");
    promptnessText.SetFillColor(0);
    promptnessText.SetTextSize(0.03);
    promptnessText.SetTextFont(62);
    promptnessText.AddText(promptness.c_str());

    TCanvas ccMc("ccMc", "ccMc", 1200, 800);
    ccMc.SetLogy(var.log_mc_);
    hmc->Draw();
    promptnessText.Draw("same");
    TPaveText disclaimer_mc(0.74, 0.78, 0.87, 0.84, "brNDC");
    disclaimer_mc.SetFillColor(0);
    disclaimer_mc.SetTextSize(0.02);
    disclaimer_mc.SetTextFont(62);
    disclaimer_mc.AddText("MC matched with reco only");
    disclaimer_mc.Draw("same");
    HistoQuantities mc_quant = EvaluateHistoQuantities(hmc);
    TPaveText mc_quant_text = ConverHistoQuantitiesToText(mc_quant, 0.74, 0.6, 0.87, 0.7);
    mc_quant_text.Draw("same");

    TCanvas ccRec("ccRec", "ccRec", 1200, 800);
    ccRec.SetLogy(var.log_rec_);
    hrec->Draw();
    promptnessText.Draw("same");
    TPaveText disclaimer_rec(0.74, 0.78, 0.87, 0.84, "brNDC");
    disclaimer_rec.SetFillColor(0);
    disclaimer_rec.SetTextSize(0.02);
    disclaimer_rec.SetTextFont(62);
    disclaimer_rec.AddText("Rec matched to MC-true only");
    disclaimer_rec.Draw("same");
    HistoQuantities rec_quant = EvaluateHistoQuantities(hrec);
    TPaveText rec_quant_text = ConverHistoQuantitiesToText(rec_quant, 0.74, 0.6, 0.87, 0.7);
    rec_quant_text.Draw("same");

    TCanvas ccRes("ccRes", "ccRes", 1200, 800);
    ccRes.SetLogy(var.log_res_);
    hres->Draw();
    promptnessText.Draw("same");
    HistoQuantities res_quant = EvaluateHistoQuantities(hres);
    TPaveText res_quant_text = ConverHistoQuantitiesToText(res_quant, 0.74, 0.6, 0.87, 0.7);
    res_quant_text.Draw("same");

    TCanvas ccCorr("ccCorr", "ccCorr", 1200, 800);
    ccCorr.SetRightMargin(0.16);
    ccCorr.SetLogz(var.log_corr_);
    hcorr->Draw("colz");
    promptnessText.Draw("same");

    TCanvas ccPull("ccPull", "ccPull", 1200, 800);
    ccPull.SetLogy(var.log_pull_);
    hpull->Draw();
    promptnessText.Draw("same");
    HistoQuantities pull_quant = EvaluateHistoQuantities(hpull);
    TPaveText pull_quant_text = ConverHistoQuantitiesToText(pull_quant, 0.74, 0.6, 0.87, 0.7);
    pull_quant_text.Draw("same");

    if(is_first_canvas) printing_bracket = "(";
    else                printing_bracket = "";
    ccMc.Print((fileOutName + "_mc.pdf" + printing_bracket).c_str(), "pdf");
    ccRec.Print((fileOutName + "_rec.pdf" + printing_bracket).c_str(), "pdf");
    ccRes.Print((fileOutName + "_res.pdf" + printing_bracket).c_str(), "pdf");
    ccCorr.Print((fileOutName + "_corr.pdf" + printing_bracket).c_str(), "pdf");
    ccPull.Print((fileOutName + "_pull.pdf" + printing_bracket).c_str(), "pdf");
    is_first_canvas = false;
  }

  printing_bracket = "]";
  TCanvas emptycanvas("emptycanvas", "", 1200, 800);
  emptycanvas.Print((fileOutName + "_mc.pdf" + printing_bracket).c_str(), "pdf");
  emptycanvas.Print((fileOutName + "_rec.pdf" + printing_bracket).c_str(), "pdf");
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
  result.overflow_ = h->GetBinContent(h->GetNbinsX()+1) / integral;
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
  text.SetTextFont(62);

  text.AddText(("underflow = " + to_string_with_precision(q.underflow_*100, 2) + "%").c_str());
  text.AddText(("overflow = " + to_string_with_precision(q.overflow_*100, 2) + "%").c_str());
  text.AddText(("#mu = " + to_string_with_precision(q.mean_, 3) + " #pm " + to_string_with_precision(q.mean_err_, 3) + " (stat.)").c_str());
  text.AddText(("#sigma = " + to_string_with_precision(q.stddev_, 3) + " #pm " + to_string_with_precision(q.stddev_err_, 3) + " (stat.)").c_str());

  return text;
}
