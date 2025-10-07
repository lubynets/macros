void ArmenterosPodolanskiCascade() {
  const std::string fileName = "/home/oleksii/alidir/working/tpc/set4/AnalysisResults.alice.data.2023.LHC23zzo.545210.apass5.2020_set4.root";

  TFile* fileIn = TFile::Open(fileName.c_str(), "read");

  TH2* hAP = fileIn->Get<TH2>("v0-selector/hCascAPplot");
  hAP->GetXaxis()->SetTitleOffset(1.75);
  hAP->GetXaxis()->SetTitle("#alpha = #frac{p_{#parallel}^{bach+(V0)} - p_{#parallel}^{V0(bach-)}}{p_{#parallel}^{bach+(V0)} + p_{#parallel}^{V0(bach-)}}");
  hAP->GetYaxis()->SetTitle("#it{q}_{T} (GeV/#it{c})");
  hAP->SetTitle("");

  const double up1 = 0.223;
  const double up2 = 0.358;
  const double up3 = 0.325;

  const double down1 = 0.197;
  const double down2 = up2;
  const double down3 = 0.25;

  const double alphaMax = 0.65;
  const double alphaMin = 0;

  const double qtMinOuterArc = 0.142;

  const double upAlphaEdge = std::sqrt(1 - qtMinOuterArc*qtMinOuterArc/up1/up1)*up3 + up2;
  const double downAlphaEdge = std::sqrt(1 - qtMinOuterArc*qtMinOuterArc/down1/down1)*down3 + down2;

  TF1* upArc = new TF1("upArc", "std::abs((std::abs(x)-[1])/[2]) < 1 && std::abs(x) >= [3] && std::abs(x) <= [4] ? [0]*TMath::Sqrt(std::abs(1 - (std::abs(x)-[1])*(std::abs(x)-[1]) / [2]/[2])) : -999", -upAlphaEdge, upAlphaEdge);
  upArc->SetParameters(up1, up2, up3, alphaMin, alphaMax);
  upArc->SetNpx(10000);

  TF1* downArc = new TF1("downArc", "std::abs((std::abs(x)-[1])/[2]) < 1 && std::abs(x) >= [3] && std::abs(x) <= [4] ? [0]*TMath::Sqrt(std::abs(1 - (std::abs(x)-[1])*(std::abs(x)-[1]) / [2]/[2])) : -999", -downAlphaEdge, downAlphaEdge);
  downArc->SetParameters(down1, down2, down3, alphaMin, alphaMax);
  downArc->SetNpx(10000);

  TLine* right = new TLine(alphaMax, downArc->Eval(alphaMax), alphaMax, upArc->Eval(alphaMax));
  TLine* left = new TLine(-alphaMax, downArc->Eval(alphaMax), -alphaMax, upArc->Eval(alphaMax));
  TLine* qtRight = new TLine(downAlphaEdge, qtMinOuterArc, upAlphaEdge, qtMinOuterArc);
  TLine* qtLeft = new TLine(-upAlphaEdge, qtMinOuterArc, -downAlphaEdge, qtMinOuterArc);

  TPaveText* xi = new TPaveText(0.80, 0.11, 0.84, 0.18);
  TPaveText* xibar = new TPaveText(-0.84, 0.11, -0.80, 0.18);
  xi->AddText("#Xi^{-}");
  xibar->AddText("#bar{#Xi^{+}}");
  xi->SetTextColor(kBlack);
  xibar->SetTextColor(kBlack);

  TPaveText* omega = new TPaveText(0.35, 0.18, 0.40, 0.28);
  TPaveText* omegabar = new TPaveText(-0.40, 0.18, -0.35, 0.28);
  omega->AddText("#Omega^{-}");
  omegabar->AddText("#bar{#Omega^{+}}");
  omega->SetTextColor(kRed);
  omegabar->SetTextColor(kRed);

  TCanvas* cc = new TCanvas("cc", "");
  gStyle->SetOptStat(0);
  cc->SetCanvasSize(1200, 800);
  cc->SetBottomMargin(0.17);
  cc->SetRightMargin(0.01);
  cc->SetTopMargin(0.02);
  cc->SetLogz();

  hAP->Draw();
  for(const auto& arc : {upArc, downArc}) {
    arc->SetLineWidth(3);
    arc->Draw("same");
  }
  for(const auto& line : {/*left, right,*/ qtRight, qtLeft}) {
    line->SetLineColor(kRed);
    line->SetLineWidth(3);
    line->Draw("same");
  }

  for(const auto& particle : {xi, xibar, omega, omegabar}) {
    particle->SetFillColor(0);
    particle->SetBorderSize(0);
    particle->SetFillStyle(0);
    particle->SetTextSize(0.08);
    particle->SetTextFont(62);
    particle->Draw("same");
  }

  cc->Print("AP.pdf", "pdf");
}
