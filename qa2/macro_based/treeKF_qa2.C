template<typename T>
std::string to_string_with_precision(const T a_value, const int n = 6);

void treeKF_qa2(const std::string& fileName) {
  TString currentMacroPath = __FILE__;
  TString directory = currentMacroPath(0, currentMacroPath.Last('/'));
  gROOT->Macro( directory + "/treeKF_qa2.style.cc" );

  TFile* fileIn = TFile::Open(fileName.c_str());
  if(fileIn == nullptr) {
    throw std::runtime_error("fileIn == nullptr");
  }

  std::string fileOutName = "treeKF_qa2.pdf";

  struct Variable {
    std::string name_;
    bool is_log_x_;
    bool is_log_y_;
  };

  std::vector<Variable> vars {
    {"Chi2prim_p",   false, false},
    {"Chi2prim_K",   false, false},
    {"Chi2prim_pi",  false, false},
    {"Chi2geo_p_pi", false, false},
    {"Chi2geo_p_K",  false, false},
    {"Chi2geo_K_pi", false, false},
    {"DCA_p_pi",     false, false},
    {"DCA_p_K",      false, false},
    {"DCA_K_pi",     false, false},
    {"Chi2geo",      false, false},
    {"Chi2topo",     false, false},
    {"LdL",          false, false},
    {"L",            false, true },
    {"T",            false, true },
  };

  struct SignalSpecies{
    std::string name_;
    Color_t color_;
  };

  std::vector<SignalSpecies> signal_species {
    {"prompt",     kRed},
//     {"nonprompt",  kBlue},
//     {"background", kGreen+2},
//     {"wrongswap",  kBlack},
//     {"data",       kBlue},
  };

  bool is_first_canvas{true};
  for(auto& var : vars) {
    std::vector<TH1D*> histos;
    for(auto& ss : signal_species) {
      histos.emplace_back(fileIn->Get<TH1D>((ss.name_ + "/h" + var.name_ + "_" + ss.name_).c_str()));
      if(histos.back() == nullptr) {
        throw std::runtime_error((ss.name_ + "/h" + var.name_ + "_" + ss.name_ + " is nullptr").c_str());
      }
    }

    TLegend* leg = new TLegend(0.76, 0.8, 0.96, 0.92);
    leg->SetBorderSize(0);
    TCanvas cc("cc", "", 1200, 800);
    cc.cd();
    cc.SetLogx(var.is_log_x_);
    cc.SetLogy(var.is_log_y_);
    double maxvalue = 0.;
    for(int iH=0; iH<histos.size(); iH++) {
      histos.at(iH)->UseCurrentStyle();
      histos.at(iH)->Scale(1./histos.at(iH)->GetEntries());
      histos.at(iH)->GetYaxis()->SetTitle("a.u.");
      histos.at(iH)->SetLineWidth(2);
      histos.at(iH)->SetLineColor(signal_species.at(iH).color_);
      if(iH == 0) histos.at(iH)->Draw("");
      else        histos.at(iH)->Draw("same");
      float underflow = histos.at(iH)->GetBinContent(0);
      float overflow = histos.at(iH)->GetBinContent(histos.at(iH)->GetNbinsX());
//       std::cout << underflow << "\t" << overflow << "\n";
      leg->AddEntry(histos.at(iH), (signal_species.at(iH).name_ + "; unfl: " + to_string_with_precision(underflow, 2) + "; ovfl: " + to_string_with_precision(overflow, 2)).c_str(), "L");
      maxvalue = std::max(maxvalue, histos.at(iH)->GetBinContent(histos.at(iH) -> GetMaximumBin()));
    }
    const double minvalue = var.is_log_y_ ? 0.0001 : 0.;
    const double maxfactor = var.is_log_y_ ? 2 : 1.1;
    histos.at(0)->GetYaxis()->SetRangeUser(minvalue, maxfactor*maxvalue);

    leg->Draw();
    if(is_first_canvas) cc.Print((fileOutName + "(").c_str(), "pdf");
    else                cc.Print(fileOutName.c_str(), "pdf");
    is_first_canvas = false;
  }

  TCanvas emptycanvas("emptycanvas", "", 1200, 800);
  emptycanvas.Print((fileOutName + "]").c_str(), "pdf");
}

template<typename T>
  std::string to_string_with_precision(const T a_value, const int n) {
  std::ostringstream out;
  out.precision(n);
  out << std::fixed << a_value;
  return out.str();
}
