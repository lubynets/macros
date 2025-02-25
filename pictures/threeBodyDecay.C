struct Arc {
  float x1_;
  float y1_;
  float x2_;
  float y2_;
  float r_;
  bool camel_;
};

TGraph* CircleFromCenter(float a, float b, float r, float phi1, float phi2);
TGraph* CircleFromArc(float x1, float y1, float x2, float y2, float r, bool camel);
TGraph* CircleFromArc(const Arc& arc);
void QuadritizeTMultiGraph(TMultiGraph* mgr);
void MarkUpCanvas(TCanvas* cc, bool is=true);
std::vector<TGraph*> CreateTubes(TGraph* gr, const Arc& arc, float tubeRadius, float compression=0.2);
TGraph* Ellips(float x, float y, float Rave, float deltaR, float alpha);

void threeBodyDecay() {
  TCanvas* cc = new TCanvas("cc", "", 1000, 1000);
  cc->SetLeftMargin(0.05);
  cc->SetBottomMargin(cc->GetLeftMargin());
  cc->SetRightMargin(0.01);
  cc->SetTopMargin(cc->GetRightMargin());

  MarkUpCanvas(cc, false);

  TMultiGraph* mgr = new TMultiGraph();

  const float Lc_x_begin = 1;
  const float Lc_y_begin = 1;
  const float Lc_x_end = 4;
  const float Lc_y_end = 2;
  const float Lc_r = 10;
  const Arc Lc_arc{Lc_x_begin, Lc_y_begin, Lc_x_end, Lc_y_end, Lc_r, true};

  const float proton_x_begin = Lc_x_end;
  const float proton_y_begin = Lc_y_end;
  const float proton_x_end = 8;
  const float proton_y_end = 2;
  const float proton_r = 10;
  const Arc proton_arc{proton_x_begin, proton_y_begin, proton_x_end, proton_y_end, proton_r, true};

  const float pion_x_begin = Lc_x_end;
  const float pion_y_begin = Lc_y_end;
  const float pion_x_end = 7;
  const float pion_y_end = -1;
  const float pion_r = 5;
  const Arc pion_arc{pion_x_begin, pion_y_begin, pion_x_end, pion_y_end, pion_r, true};

  const float kaon_x_begin = Lc_x_end;
  const float kaon_y_begin = Lc_y_end;
  const float kaon_x_end = 6;
  const float kaon_y_end = 4;
  const float kaon_r = 7;
  const Arc kaon_arc{kaon_x_begin, kaon_y_begin, kaon_x_end, kaon_y_end, kaon_r, false};

  auto Lc_sim = CircleFromArc(Lc_arc);
  Lc_sim->SetLineColor(kOrange-7);
  Lc_sim->SetLineWidth(3);
  Lc_sim->SetLineStyle(1);
  mgr->Add(Lc_sim, "L");

  auto Lc_tubes = CreateTubes(Lc_sim, Lc_arc, 0.1);
  for(auto& tube : Lc_tubes) {
    tube->SetLineColor(Lc_sim->GetLineColor());
    tube->SetLineWidth(1);
    tube->SetLineStyle(1);
    mgr->Add(tube, "L");
  }

  auto proton_sim = CircleFromArc(proton_arc);
  proton_sim->SetLineColor(kRed);
  proton_sim->SetLineWidth(3);
  proton_sim->SetLineStyle(1);
  mgr->Add(proton_sim, "L");

  auto proton_tubes = CreateTubes(proton_sim, proton_arc, 0.1);
  for(auto& tube : proton_tubes) {
    tube->SetLineColor(proton_sim->GetLineColor());
    tube->SetLineWidth(1);
    tube->SetLineStyle(1);
    mgr->Add(tube, "L");
  }

  auto pion_sim = CircleFromArc(pion_arc);
  pion_sim->SetLineColor(kBlue);
  pion_sim->SetLineWidth(3);
  pion_sim->SetLineStyle(1);
  mgr->Add(pion_sim, "L");

  auto pion_tubes = CreateTubes(pion_sim, pion_arc, 0.1);
  for(auto& tube : pion_tubes) {
    tube->SetLineColor(pion_sim->GetLineColor());
    tube->SetLineWidth(1);
    tube->SetLineStyle(1);
    mgr->Add(tube, "L");
  }

  auto kaon_sim = CircleFromArc(kaon_arc);
  kaon_sim->SetLineColor(kGreen+2);
  kaon_sim->SetLineWidth(3);
  kaon_sim->SetLineStyle(1);
  mgr->Add(kaon_sim, "L");

  auto kaon_tubes = CreateTubes(kaon_sim, kaon_arc, 0.1);
  for(auto& tube : kaon_tubes) {
    tube->SetLineColor(kaon_sim->GetLineColor());
    tube->SetLineWidth(1);
    tube->SetLineStyle(1);
    mgr->Add(tube, "L");
  }

  TGraph* verteces = new TGraph();
  verteces->AddPoint(Lc_x_begin, Lc_y_begin);
  verteces->AddPoint(Lc_x_end, Lc_y_end);
  verteces->SetMarkerStyle(kFullCircle);
  verteces->SetMarkerSize(1);
  mgr->Add(verteces, "P");

  QuadritizeTMultiGraph(mgr);
  mgr->Draw("A");
}

std::vector<TGraph*> CreateTubes(TGraph* gr, const Arc& arc, float tubeRadius, float compression) {
  const float xFirst = gr->GetPointX(0);
  const float yFirst = gr->GetPointY(0);
  const float xSecond = gr->GetPointX(1);
  const float ySecond = gr->GetPointY(1);
  const float xLast = gr->GetPointX(gr->GetN()-1);
  const float yLast = gr->GetPointY(gr->GetN()-1);
  const float xPreLast = gr->GetPointX(gr->GetN()-2);
  const float yPreLast = gr->GetPointY(gr->GetN()-2);

  const float phi1 = TMath::ATan2(ySecond-yFirst, xSecond-xFirst);
  const float phi2 = TMath::ATan2(yLast-yPreLast, xLast-xPreLast);

  const float phi1Rot = phi1 + TMath::Pi()/2;
  const float phi2Rot = phi2 + TMath::Pi()/2;

  const float x1Up = xFirst + tubeRadius * TMath::Cos(phi1Rot);
  const float x1Down = xFirst - tubeRadius * TMath::Cos(phi1Rot);
  const float y1Up = yFirst + tubeRadius * TMath::Sin(phi1Rot);
  const float y1Down = yFirst - tubeRadius * TMath::Sin(phi1Rot);

  const float x2Up = xLast + tubeRadius * TMath::Cos(phi2Rot);
  const float x2Down = xLast - tubeRadius * TMath::Cos(phi2Rot);
  const float y2Up = yLast + tubeRadius * TMath::Sin(phi2Rot);
  const float y2Down = yLast - tubeRadius * TMath::Sin(phi2Rot);

  const double q = arc.camel_ ? 1 : -1;

  TGraph* tubeUp = CircleFromArc(x1Up, y1Up, x2Up, y2Up, arc.r_+q*tubeRadius, arc.camel_);
  TGraph* tubeDown = CircleFromArc(x1Down, y1Down, x2Down, y2Down, arc.r_-q*tubeRadius, arc.camel_);

  const float Rave = tubeRadius*(1+compression) / 2;
  const float deltaR = tubeRadius*(1-compression) / 2;

  TGraph* ellips1 = Ellips(xFirst, yFirst, Rave, deltaR, phi1Rot);
  TGraph* ellips2 = Ellips(xLast, yLast, Rave, deltaR, phi2Rot);

  return {tubeUp, tubeDown, ellips1, ellips2};
}

TGraph* CircleFromCenter(float a, float b, float r, float phi1, float phi2) {
  const int N = 1000;
  TGraph* gr = new TGraph();
  for (int i=0; i<N; i++) {
    const float phi = phi1 + (phi2-phi1)*i/N;
    const float x = a + r*TMath::Cos(phi);
    const float y = b + r*TMath::Sin(phi);
    gr->AddPoint(x,y);
  }

  return gr;
}

TGraph* CircleFromArc(const Arc& arc) {
  return CircleFromArc(arc.x1_, arc.y1_, arc.x2_, arc.y2_, arc.r_, arc.camel_);
}

TGraph* CircleFromArc(float x1, float y1, float x2, float y2, float r, bool camel) {
  const float d2 = (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1);
  const float d = std::sqrt(d2);
  const float D2 = r*r - d2/4;
  if(D2 < 0) throw std::runtime_error("CircleFromArc(): circle points and its radius are incompatible");
  const float D = std::sqrt(D2);

  const float q = camel ? 1 : -1;

  const float xc = (x1+x2)/2 + q * D * (y2-y1)/d;
  const float yc = (y1+y2)/2 - q * D * (x2-x1)/d;

  const float phi1 = TMath::ATan2(y1-yc, x1-xc);
  const float phi2 = TMath::ATan2(y2-yc, x2-xc);

  return CircleFromCenter(xc, yc, r, phi1, phi2);
}

void QuadritizeTMultiGraph(TMultiGraph* mgr) {
  const double xmin = mgr->GetXaxis()->GetXmin();
  const double xmax = mgr->GetXaxis()->GetXmax();
  const double ymin = mgr->GetYaxis()->GetXmin();
  const double ymax = mgr->GetYaxis()->GetXmax();

  const double min = std::min(xmin, ymin);
  const double max = std::max(xmax, ymax);

  mgr->GetXaxis()->SetRangeUser(min, max);
  mgr->GetYaxis()->SetRangeUser(min, max);
}

TGraph* Ellips(float xc, float yc, float Rave, float deltaR, float alpha) {
  const float a = Rave + deltaR;
  const float b = Rave - deltaR;
  if(a < b) throw std::runtime_error("Ellips(): a < b");
  const float e = std::sqrt(a*a - b*b) / a;
  const float f = a*e;
  const float xref = xc - f*TMath::Cos(alpha);
  const float yref = yc - f*TMath::Sin(alpha);

  const int N = 1000;
  TGraph* gr = new TGraph();
  for (int i=0; i<N; i++) {
    const float phi = 2*TMath::Pi() * i / N;
    const float rho = - a * (1-e*e) / (1 + e*TMath::Cos(-alpha + phi));
    const float x = xref + rho * TMath::Cos(phi);
    const float y = yref + rho * TMath::Sin(phi);
    gr->AddPoint(x, y);
  }

  return gr;
}

void MarkUpCanvas(TCanvas* cc, bool is) {
  if(is) {
    cc->SetGrid();
  } else {
    gStyle->SetTickLength(0.0,"xy");
  }
}
