void pt_gen_builder(const std::string fileListName) {
  const std::string treeName = "aTree";
  const std::string ptGenFieldName = "fGen_Pt";

  AnalysisTree::Chain* treeIn = new AnalysisTree::Chain((std::vector<std::string>){fileListName}, (std::vector<std::string>){treeName});

  AnalysisTree::GenericDetector* genParticles{nullptr};

  const int ptGenFieldId = treeIn->GetConfiguration()->GetBranchConfig("Generated").GetFieldId("fGen_Pt");

  treeIn->SetBranchAddress("Generated.", &genParticles);

  TH1D* histoPtGen = new TH1D("histoPtGen", "", 2000, 0, 20);
  histoPtGen->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  histoPtGen->GetYaxis()->SetTitle("Entries");

  const int nEntries = treeIn->GetEntries();
  for(int iEntry=0; iEntry<nEntries; ++iEntry) {
    treeIn->GetEntry(iEntry);
    for(const auto& genParticle : *(genParticles->GetChannels()) ) {
      histoPtGen->Fill(genParticle.GetField<float>(ptGenFieldId));
    }
  }

  histoPtGen->SaveAs("ptGen.root");
}
