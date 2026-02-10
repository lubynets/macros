void plot_histograms(const std::string& fileName, const std::string& dirName) {
  gStyle->SetOptStat("nemruo");
  TFile* f = TFile::Open(fileName.c_str(), "read");
  if (!f || f->IsZombie()) {
    std::cerr << "Cannot open file\n";
    return;
  }

  // Access directory
  TDirectory* dir = f->GetDirectory(dirName.c_str());
  if (!dir) {
    std::cerr << "Directory " << dirName << " not found\n";
    return;
  }

  // Prepare canvas and PDF
  TCanvas canvas("c", "Histograms", 1200, 800);
  canvas.SetLogz();
  std::string pdfName;
  if(dirName.empty()) {
    pdfName = "histograms.pdf";
  } else {
    pdfName = dirName + ".pdf";
    std::replace(pdfName.begin(), pdfName.end(), '/', '_');
  }
  canvas.Print(Form("%s[", pdfName.c_str()));  // Start multi-page PDF

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
    canvas.Print(pdfName.c_str());  // Add page to PDF
  }

  canvas.Print(Form("%s]", pdfName.c_str()));  // End multi-page PDF

  f->Close();
}
