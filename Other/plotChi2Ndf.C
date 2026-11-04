void plotChi2Ndf(int ndf) {
  TF1* f = new TF1("chi2", [ndf](double* x, double*) { return ROOT::Math::chisquared_pdf(x[0], ndf); }, 0., 20, 0);

  f->SetTitle(Form("#chi^{2} distribution, ndf = %d;#chi^{2};PDF", ndf));
  f->Draw();
}
