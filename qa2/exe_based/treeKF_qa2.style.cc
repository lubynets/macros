//
// Created by mikhail on 3/5/21.
//
{
gStyle->SetCanvasPreferGL(true);
gStyle->SetPadLeftMargin(0.14);
gStyle->SetPadRightMargin(0.02);
gStyle->SetPadBottomMargin(0.12);
gStyle->SetPadTopMargin(0.07);
gStyle->SetLegendBorderSize(0);
gStyle->SetFrameLineWidth(4);
gStyle->SetMarkerSize(2);
gStyle->SetLineWidth(3);
gStyle->SetHistLineWidth(3);

gStyle->SetTitleSize(0.05, "X");
gStyle->SetTitleSize(0.05, "Y");
gStyle->SetTitleSize(0.01, "Z");

gStyle->SetLabelSize(0.05, "X");
gStyle->SetLabelSize(0.05, "Y");
gStyle->SetLabelSize(0.05, "Z");

gStyle->SetLabelOffset(0.003, "X");
gStyle->SetLabelOffset(0.01, "Y");
gStyle->SetLabelOffset(0.05, "Z");

gStyle->SetTitleOffset(1.0, "X");
gStyle->SetTitleOffset(1.46, "Y");
gStyle->SetTitleOffset(0.05, "Z");

gStyle->SetNdivisions(206, "xyz");

gStyle->SetOptStat(0);
}
