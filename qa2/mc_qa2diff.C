template<typename T>
std::string to_string_with_precision(const T a_value, const int n = 6);

std::vector<std::pair<std::string, std::string>> FindCuts(const TFile* fileIn, std::string name_start);

bool stofCompare(std::pair<std::string, std::string> a, std::pair<std::string, std::string> b);

struct HistoQuantities {
  float underflow_{-999.f};
  float overflow_{-999.f};
  float mean_{-999.f};
  float mean_err_{-999.f};
  float stddev_{-999.f};
  float stddev_err_{-999.f};
};

HistoQuantities EvaluateHistoQuantities(const TH1* h);
TPaveText ConvertHistoQuantitiesToText(const HistoQuantities& q, float x1, float y1, float x2, float y2);

void mc_qa2diff(const std::string& fileName, int prompt_or_nonprompt=1) {
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
//  name    cutname  title          unit
    {"P",   "psim",  "p^{mc}",     "GeV/c", false, false},
    {"Pt",  "pTsim", "p_{T}^{mc}", "GeV/c", false, false},
    {"Xsv", "psim",  "p^{mc}",     "GeV/c", false, false},
    {"Ysv", "psim",  "p^{mc}",     "GeV/c", false, false},
    {"Zsv", "psim",  "p^{mc}",     "GeV/c", false, false},
    {"L",   "lsim",  "L^{mc}",     "cm",    true,  true},
    {"T",   "tsim",  "T^{mc}",     "ps",    true,  true},
  };

  for(auto& var : vars) {
    bool is_first_canvas{true};
    std::string printing_bracket = "(";

    auto cutsVar = FindCuts(fileIn, "Candidates_Simulated_"  + promptness + "_" + var.cut_name_);

    TGraphMultiErrors grResMu(cutsVar.size(), 2);
    grResMu.SetTitle("");
    grResMu.GetXaxis()->SetTitle((var.cut_title_ + ", " + var.cut_unit_).c_str());
    grResMu.SetMarkerStyle(kFullSquare);
    grResMu.SetMarkerSize(1.6);
    grResMu.SetMarkerColor(kBlue);
    grResMu.SetLineWidth(3);
    grResMu.SetLineColor(kBlue);

    int iPoint{0};
    for(auto& cV : cutsVar) {
      TH1D* hRes = fileIn->Get<TH1D>(("Candidates_Simulated_"  + promptness + "_"  + var.cut_name_ + "_"  + cV.first + "_"  + cV.second + "/res_"  + var.name_ + "_"  + var.cut_name_ + "_"  + cV.first + "_"  + cV.second).c_str());
      TH1D* hPull = fileIn->Get<TH1D>(("Candidates_Simulated_"  + promptness + "_"  + var.cut_name_ + "_"  + cV.first + "_"  + cV.second + "/res_"  + var.name_ + "_"  + var.cut_name_ + "_"  + cV.first + "_"  + cV.second).c_str());
      if(hRes == nullptr || hPull == nullptr) {
        throw std::runtime_error("hhRes == nullptr || hPull == nullptr for " + var.name_ + "; " + var.cut_name_ + "; " + cV.first + "; "  + cV.second);
      }

      hRes->UseCurrentStyle();
      hPull->UseCurrentStyle();

      if(is_first_canvas) {
        grResMu.GetYaxis()->SetTitle(("#mu_{" + static_cast<std::string>(hRes->GetXaxis()->GetTitle()) + "}").c_str());
      }

      TPaveText generalText(0.74, 0.82, 0.87, 0.90, "brNDC");
      generalText.SetFillColor(0);
      generalText.SetTextSize(0.03);
      generalText.SetTextFont(62);
      generalText.AddText(promptness.c_str());
      generalText.AddText((cV.first + " < " + var.cut_title_ + " < " + cV.second + " " + var.cut_unit_).c_str());

      TCanvas ccRes("ccRes", "ccRes", 1200, 800);
      ccRes.SetLogy(var.log_res_);
      hRes->Draw();
      generalText.Draw("same");
      HistoQuantities res_quant = EvaluateHistoQuantities(hRes);
      TPaveText res_quant_text = ConvertHistoQuantitiesToText(res_quant, 0.74, 0.6, 0.87, 0.7);
      res_quant_text.Draw("same");

      const float X = (atof(cV.first.c_str()) + atof(cV.second.c_str())) / 2;
      grResMu.SetPoint(iPoint, X, res_quant.mean_);
      const float eX = /*(atof(cV.second.c_str()) - atof(cV.first.c_str())) / 1*/0;
      grResMu.SetPointEX(iPoint, eX, eX);
      grResMu.SetPointEY(iPoint, 0, res_quant.mean_err_, res_quant.mean_err_);
      grResMu.SetPointEY(iPoint, 1, res_quant.stddev_, res_quant.stddev_);

      TCanvas ccPull("ccPull", "ccPull", 1200, 800);
      ccPull.SetLogy(var.log_pull_);
      hPull->Draw();
      generalText.Draw("same");
      HistoQuantities pull_quant = EvaluateHistoQuantities(hPull);
      TPaveText pull_quant_text = ConvertHistoQuantitiesToText(pull_quant, 0.74, 0.6, 0.87, 0.7);
      pull_quant_text.Draw("same");

      if(is_first_canvas) printing_bracket = "(";
      else                printing_bracket = "";
      ccRes.Print((fileOutName + "_" + var.name_ + "_res.pdf" + printing_bracket).c_str(), "pdf");
      ccPull.Print((fileOutName + "_" + var.name_ + "_pull.pdf" + printing_bracket).c_str(), "pdf");
      is_first_canvas = false;
      ++iPoint;
    } // cutsVar

    TCanvas ccResMu("ccResMu", "ccResMu", 1200, 800);
    grResMu.Draw("AP; []");
    ccResMu.Print((fileOutName + "_" + var.name_ + "_res.pdf" + printing_bracket).c_str(), "pdf");

    printing_bracket = "]";
    TCanvas emptycanvas("emptycanvas", "", 1200, 800);
    emptycanvas.Print((fileOutName + "_" + var.name_ + "_res.pdf" + printing_bracket).c_str(), "pdf");
    emptycanvas.Print((fileOutName + "_" + var.name_ + "_pull.pdf" + printing_bracket).c_str(), "pdf");
  } // vars

  fileIn->Close();
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

TPaveText ConvertHistoQuantitiesToText(const HistoQuantities& q, float x1, float y1, float x2, float y2) {
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

  return result;
}

bool stofCompare(std::pair<std::string, std::string> a, std::pair<std::string, std::string> b) {
  return atof(a.first.c_str()) < atof(b.first.c_str());
}
