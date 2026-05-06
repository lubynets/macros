void compareEfficiencyTwoMcs() {
  const bool isPlotRatio{false};
  const bool isDrawNoW{false};

  const std::vector<std::string> npScores{"0.20", "0.40", "0.70"};
  const std::vector<Color_t> palette{kRed, kGreen+2, kBlue, kOrange+2, kMagenta};

  TFile* fileOld = TFile::Open("/home/oleksii/alidir/working/systematics/ctMc/HF_LHC24h1b_All/667716/effWnew/efficiency_summary.root");
  TFile* fileNew = TFile::Open("/home/oleksii/alidir/working/systematics/ctMc/HF_LHC24h1c_All/669133/effWnew/efficiency_summary.root");

  TCanvas* cPrompt = new TCanvas("cPrompt", "", 1200, 800);
  TCanvas* cNonPrompt = new TCanvas("cNonPrompt", "", 1200, 800);

  TMultiGraph* mgPrompt = new TMultiGraph();
  TMultiGraph* mgNonPrompt = new TMultiGraph();

  const std::string yAxisTitle = isPlotRatio ? "#varepsilon [h1c] / #varepsilon [h1b]" : "#varepsilon";

  mgPrompt->SetTitle(("Prompt;#it{t} (ps);" + yAxisTitle).c_str());
  mgNonPrompt->SetTitle(("NonPrompt;#it{t} (ps);" + yAxisTitle).c_str());

  TLegend* leg = new TLegend(0.75, 0.75, 0.99, 0.99);

  const Style_t oldStyle = kFullSquare;
  const Style_t newStyle = kFullCircle;

  if(!isPlotRatio) {
    auto* entry = leg->AddEntry("", "HF_LHC24h1b_All", "P");
    entry->SetMarkerSize(2);
    entry->SetMarkerColor(kBlack);
    entry->SetMarkerStyle(oldStyle);
    entry = leg->AddEntry("", "HF_LHC24h1c_All", "P");
    entry->SetMarkerSize(2);
    entry->SetMarkerColor(kBlack);
    entry->SetMarkerStyle(newStyle);
  }

  for(int iNp=0, nNps=npScores.size(); iNp<nNps; ++iNp) {
    TH1* hOldPromptNoW = fileOld->Get<TH1>(("effs/prompt/pT_3_20/eff_NPgt" + npScores.at(iNp)).c_str());
    TH1* hOldPromptW = fileOld->Get<TH1>(("effs/prompt/pT_3_20/eff_NPgt" + npScores.at(iNp) + "_W").c_str());
    TH1* hOldNonPromptNoW = fileOld->Get<TH1>(("effs/nonprompt/pT_3_20/eff_NPgt" + npScores.at(iNp)).c_str());
    TH1* hOldNonPromptW = fileOld->Get<TH1>(("effs/nonprompt/pT_3_20/eff_NPgt" + npScores.at(iNp) + "_W").c_str());
    TH1* hNewPromptNoW = fileNew->Get<TH1>(("effs/prompt/pT_3_20/eff_NPgt" + npScores.at(iNp)).c_str());
    TH1* hNewPromptW = fileNew->Get<TH1>(("effs/prompt/pT_3_20/eff_NPgt" + npScores.at(iNp) + "_W").c_str());
    TH1* hNewNonPromptNoW = fileNew->Get<TH1>(("effs/nonprompt/pT_3_20/eff_NPgt" + npScores.at(iNp)).c_str());
    TH1* hNewNonPromptW = fileNew->Get<TH1>(("effs/nonprompt/pT_3_20/eff_NPgt" + npScores.at(iNp) + "_W").c_str());

    if(isPlotRatio) {
      hNewPromptNoW->Divide(hOldPromptNoW);
      hNewPromptW->Divide(hOldPromptW);
      hNewNonPromptNoW->Divide(hOldNonPromptNoW);
      hNewNonPromptW->Divide(hOldNonPromptW);
    }

    TGraph* grOldPromptNoW = new TGraph(hOldPromptNoW);
    TGraph* grOldPromptW = new TGraph(hOldPromptW);
    TGraph* grOldNonPromptNoW = new TGraph(hOldNonPromptNoW);
    TGraph* grOldNonPromptW = new TGraph(hOldNonPromptW);
    TGraph* grNewPromptNoW = new TGraph(hNewPromptNoW);
    TGraph* grNewPromptW = new TGraph(hNewPromptW);
    TGraph* grNewNonPromptNoW = new TGraph(hNewNonPromptNoW);
    TGraph* grNewNonPromptW = new TGraph(hNewNonPromptW);

    for(const auto& gr : {grOldPromptNoW, grOldPromptW, grOldNonPromptNoW, grOldNonPromptW, grNewPromptNoW, grNewPromptW, grNewNonPromptNoW, grNewNonPromptW}) {
      gr->SetMarkerColor(palette.at(iNp));
      gr->SetMarkerSize(2);
    }

    for(const auto& gr : {grOldPromptNoW, grOldNonPromptNoW}) gr->SetMarkerStyle(kOpenSquare);
    for(const auto& gr : {grOldPromptW, grOldNonPromptW}) gr->SetMarkerStyle(oldStyle);
    for(const auto& gr : {grNewPromptNoW, grNewNonPromptNoW}) gr->SetMarkerStyle(kOpenCircle);
    for(const auto& gr : {grNewPromptW, grNewNonPromptW}) gr->SetMarkerStyle(newStyle);

    mgPrompt->Add(grNewPromptW, "P");
    mgNonPrompt->Add(grNewNonPromptW, "P");
    if(isDrawNoW) {
      mgPrompt->Add(grNewPromptNoW, "P");
      mgNonPrompt->Add(grNewNonPromptNoW, "P");
    }

    if(!isPlotRatio) {
      mgPrompt->Add(grOldPromptW, "P");
      mgNonPrompt->Add(grOldNonPromptW, "P");
      if(isDrawNoW) {
        mgPrompt->Add(grOldPromptNoW, "P");
        mgNonPrompt->Add(grOldNonPromptNoW, "P");
      }
    }

    leg->AddEntry(grNewPromptW, ("bdt(np) > " + npScores.at(iNp)).c_str(), "P");
  }

  cPrompt->cd();
  mgPrompt->Draw("AP");
  leg->Draw("same");
  cPrompt->Print("effs.pdf(", "pdf");

  cNonPrompt->cd();
  mgNonPrompt->Draw("AP");
  leg->Draw("same");
  cNonPrompt->Print("effs.pdf)", "pdf");
}
