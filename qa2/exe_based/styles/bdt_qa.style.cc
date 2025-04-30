{
gStyle->SetCanvasPreferGL(true);
gStyle->SetPadLeftMargin(0.12);
gStyle->SetPadRightMargin(0.17);
gStyle->SetPadBottomMargin(0.09);
gStyle->SetPadTopMargin(0.07);
gStyle->SetLegendBorderSize(0);
gStyle->SetFrameLineWidth(4);
gStyle->SetMarkerSize(2);
gStyle->SetLineWidth(3);
gStyle->SetHistLineWidth(3);

gStyle->SetTitleSize(0.04, "X");
gStyle->SetTitleSize(0.04, "Y");
gStyle->SetTitleSize(0.04, "Z");

gStyle->SetLabelSize(0.04, "X");
gStyle->SetLabelSize(0.04, "Y");
gStyle->SetLabelSize(0.04, "Z");

gStyle->SetLabelOffset(0.003, "X");
gStyle->SetLabelOffset(0.01, "Y");
gStyle->SetLabelOffset(0.01, "Z");

gStyle->SetTitleOffset(1.0, "X");
gStyle->SetTitleOffset(1.46, "Y");
gStyle->SetTitleOffset(1.05, "Z");

gStyle->SetNdivisions(206, "xyz");

gStyle->SetOptStat(0);
}