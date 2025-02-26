#include "DecayDrawer.h"

void threeBodyDecay() {
  TCanvas* cc = new TCanvas("cc", "", 1200, 800);
  cc->SetLeftMargin(0.05);
  cc->SetBottomMargin(cc->GetLeftMargin());
  cc->SetRightMargin(0.01);
  cc->SetTopMargin(cc->GetRightMargin());

  MarkUpCanvas(cc, false);

  TMultiGraph* mgr = new TMultiGraph();

  // primary and secondary vertex
  const float pvX = 1;
  const float pvY = 1;
  const float svX = 4;
  const float svY = 2;

  // Lambda c
  const float Lc_begin_shift_X = -0.1;
  const float Lc_begin_shift_Y = 0.08;
  const float Lc_end_shift_X = -0.1;
  const float Lc_end_shift_Y = 0.06;

  const float Lc_x_begin = pvX + Lc_begin_shift_X;
  const float Lc_y_begin = pvY + Lc_begin_shift_Y;
  const float Lc_x_end = svX + Lc_end_shift_X;
  const float Lc_y_end = svY + Lc_end_shift_Y;
  const float Lc_r = 10;
  const Arc Lc_arc{Lc_x_begin, Lc_y_begin, Lc_x_end, Lc_y_end, Lc_r, true};

  // proton
  const float proton_begin_shift_X = 0.12;
  const float proton_begin_shift_Y = 0.04;

  const float proton_x_begin = svX + proton_begin_shift_X;
  const float proton_y_begin = svY + proton_begin_shift_Y;
  const float proton_x_end = 8;
  const float proton_y_end = 2;
  const float proton_r = 10;
  const Arc proton_arc{proton_x_begin, proton_y_begin, proton_x_end, proton_y_end, proton_r, true};

  // pion
  const float pion_begin_shift_X = -0.09;
  const float pion_begin_shift_Y = -0.11;

  const float pion_x_begin = svX + pion_begin_shift_X;
  const float pion_y_begin = svY + pion_begin_shift_Y;
  const float pion_x_end = 7;
  const float pion_y_end = -1;
  const float pion_r = 5;
  const Arc pion_arc{pion_x_begin, pion_y_begin, pion_x_end, pion_y_end, pion_r, true};

  // kaon
  const float kaon_begin_shift_X = -0.04;
  const float kaon_begin_shift_Y = 0.1;

  const float kaon_x_begin = svX + kaon_begin_shift_X;
  const float kaon_y_begin = svY + kaon_begin_shift_Y;
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
  verteces->AddPoint(pvX, pvY);
  verteces->AddPoint(svX, svY);
  verteces->SetMarkerStyle(kFullCircle);
  verteces->SetMarkerSize(1);
  mgr->Add(verteces, "P");

  mgr->Draw("A");
}
