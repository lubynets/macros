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
TGraph* Ellipse(float x, float y, float Rave, float deltaR, float alpha);

// ===========================================================================================================

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

  TGraph* ellipse1 = Ellipse(xFirst, yFirst, Rave, deltaR, phi1Rot);
  TGraph* ellipse2 = Ellipse(xLast, yLast, Rave, deltaR, phi2Rot);

  return {tubeUp, tubeDown, ellipse1, ellipse2};
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

TGraph* Ellipse(float xc, float yc, float Rave, float deltaR, float alpha) {
  const float a = Rave + deltaR;
  const float b = Rave - deltaR;
  if(a < b) throw std::runtime_error("Ellipse(): a < b");
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
