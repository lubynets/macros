void plot_histograms(const std::string& fileName, const std::string& dirName) {
  TFile* f = TFile::Open(fileName.c_str(), "read");
  if (!f || f->IsZombie()) {
    std::cerr << "Cannot open file\n";
    return;
  }

  // Access directory
  TDirectory* dir = f->GetDirectory(dirName.c_str());
  if (!dir) {
    std::cerr << "Directory 'v0-selector' not found\n";
    return;
  }

  // Prepare canvas and PDF
  TCanvas canvas("c", "Histograms", 1200, 800);
  canvas.SetLogz();
  const char* pdfName = "histograms.pdf";
  canvas.Print(Form("%s[", pdfName));  // Start multi-page PDF

  // Iterate over objects in directory
  TIter keyIter(dir->GetListOfKeys());
  TKey* key;
  while ((key = (TKey*)keyIter())) {
    TClass* cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1")) continue;  // Skip if not a histogram

    TH1* h = (TH1*)key->ReadObj();
    if (!h) continue;

    canvas.Clear();
    h->Draw();
    canvas.Print(pdfName);  // Add page to PDF
  }

  canvas.Print(Form("%s]", pdfName));  // End multi-page PDF

  f->Close();
}
