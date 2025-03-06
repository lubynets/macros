void AddExperiment(TGraphErrors* gr, float y, std::vector<float> t);
void AddText(const std::string& text, float x, float y);

void lifetimes() {
  TString currentMacroPath = __FILE__;
  TString directory = currentMacroPath(0, currentMacroPath.Last('/'));
  gROOT->Macro( directory + "/lifetimes.style.cc" );

  struct Experiment {
    std::string name_;
    std::string note_;
    float y_value_;
  };

  std::vector<Experiment> experiments {
    {"PDG",      "(2018)",    4},
    {"LHCb",     "(2018/19)", 3},
    {"LHCb",     "(2021)",    2},
    {"Belle II", "(2023)",    1}
  };

  struct Particle {
    std::string name_;
    std::string greek_name_;
    Color_t color_;
    std::vector<float> pdg_;
    std::vector<float> lhcb_cl_;
    std::vector<float> lhcb_prompt_;
    std::vector<float> belleii_;
  };

  std::vector<Particle> particles {//             pdg         lhcb sl             lhcb prompt         belle ii
    {"Lambda_c", "#Lambda^{+}_{c}", kViolet,   {200, 6, 0},  {202.1, 1.7, 0.9}, {},                 {203.2, 0.89, 0.77}},
    {"Xi_c0",    "#Xi^{0}_{c}",     kOrange+2, {112, 13, 0}, {153.4, 2.4, 0.7}, {148, 2.3, 2.2},    {}},
    {"Xi_c+",    "#Xi^{+}_{c}",     kBlue,     {442, 26, 0}, {454, 5, 2},       {},                 {}},
    {"Omega_c0", "#Omega^{0}_{c}",  kRed,      {69, 12, 0},  {268, 24, 10},     {276.5, 13.4, 4.5}, {243, 48, 11}},
  };

  TCanvas* cc = new TCanvas("cc", "", 1200, 800);
  for(auto& pa : particles) {
    TGraphErrors* gr = new TGraphErrors();
    AddExperiment(gr, experiments.at(0).y_value_, pa.pdg_);
    AddExperiment(gr, experiments.at(1).y_value_, pa.lhcb_cl_);
    AddExperiment(gr, experiments.at(2).y_value_, pa.lhcb_prompt_);
    if(pa.name_ != "Lambda_c") AddExperiment(gr, experiments.at(3).y_value_, pa.belleii_);
    else                       AddExperiment(gr, experiments.at(3).y_value_+0.2, pa.belleii_);
    gr->SetMarkerColor(pa.color_);
    gr->SetLineColor(pa.color_);
    gr->SetMarkerSize(2);
    gr->SetMarkerStyle(kFullSquare);
    if(pa.name_ == "Lambda_c") {
      gr->GetXaxis()->SetTitle("lifetime (fs)");
      gr->GetXaxis()->SetLimits(0, 500);
      gr->GetYaxis()->Set(40, 0, 5);
      gr->GetYaxis()->SetTickSize(0);
      for(int iExp=0; iExp<experiments.size(); iExp++) {
        gr->GetYaxis()->SetBinLabel(gr->GetYaxis()->FindBin(experiments.at(iExp).y_value_+0.125), experiments.at(iExp).name_.c_str());
        gr->GetYaxis()->SetBinLabel(gr->GetYaxis()->FindBin(experiments.at(iExp).y_value_-0.125), experiments.at(iExp).note_.c_str());
      }
      gr->GetYaxis()->SetRangeUser(0, 5);
      gr->Draw("AP");
    } else {
      gr->Draw("P");
    }

    AddText(particles.at(0).greek_name_, particles.at(0).pdg_.at(0), experiments.at(0).y_value_ + 0.3);
    AddText(particles.at(0).greek_name_, particles.at(0).lhcb_cl_.at(0), experiments.at(1).y_value_ + 0.3);
//     AddText(particles.at(0).greek_name_, particles.at(0).lhcb_prompt_.at(0), experiments.at(2).y_value_ + 0.3);
    AddText(particles.at(0).greek_name_, particles.at(0).belleii_.at(0), experiments.at(3).y_value_ + 0.3 + 0.2);

    AddText(particles.at(1).greek_name_, particles.at(1).pdg_.at(0), experiments.at(0).y_value_ + 0.3);
    AddText(particles.at(1).greek_name_, particles.at(1).lhcb_cl_.at(0), experiments.at(1).y_value_ + 0.3);
    AddText(particles.at(1).greek_name_, particles.at(1).lhcb_prompt_.at(0), experiments.at(2).y_value_ + 0.3);
//     AddText(particles.at(1).greek_name_, particles.at(1).belleii_.at(0), experiments.at(3).y_value_ + 0.3);

    AddText(particles.at(2).greek_name_, particles.at(2).pdg_.at(0), experiments.at(0).y_value_ + 0.3);
    AddText(particles.at(2).greek_name_, particles.at(2).lhcb_cl_.at(0), experiments.at(1).y_value_ + 0.3);
//     AddText(particles.at(2).greek_name_, particles.at(2).lhcb_prompt_.at(0), experiments.at(2).y_value_ + 0.3);
//     AddText(particles.at(2).greek_name_, particles.at(2).belleii_.at(0), experiments.at(3).y_value_ + 0.3);

    AddText(particles.at(3).greek_name_, particles.at(3).pdg_.at(0), experiments.at(0).y_value_ + 0.3);
    AddText(particles.at(3).greek_name_, particles.at(3).lhcb_cl_.at(0), experiments.at(1).y_value_ + 0.3);
    AddText(particles.at(3).greek_name_, particles.at(3).lhcb_prompt_.at(0), experiments.at(2).y_value_ + 0.3);
    AddText(particles.at(3).greek_name_, particles.at(3).belleii_.at(0), experiments.at(3).y_value_ + 0.3);
  }

  cc->Print("lifetimes.pdf", "pdf");

}

void AddText(const std::string& text, float x, float y) {
  const float w = 0.05;
  const float h = 0.05;
  TPaveText *pt2 = new TPaveText(x-w/2, y-h/2, x+w/2, y+h/2);
  pt2->SetFillColor(0);
  pt2->SetBorderSize(0);
  pt2->SetFillStyle(0);
  pt2->AddText(text.c_str());
  pt2->SetTextSize(0.05);
  pt2->SetTextFont(62);
  pt2->Draw("same");
}

void AddExperiment(TGraphErrors* gr, float y, std::vector<float> t) {
  if(t.size() == 0) return;
  const float x = t.at(0);
  const float ex = std::sqrt(t.at(1)*t.at(1) + t.at(2)*t.at(2));
  const int iP = gr->GetN();
  gr->SetPoint(iP, x, y);
  gr->SetPointError(iP, ex, 0);
}
