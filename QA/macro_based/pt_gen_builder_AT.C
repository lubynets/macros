void pt_gen_builder_AT(const std::string fileListName) {
  const std::string treeName = "aTree";
  const std::string ptGenFieldName = "fGen_Pt";

  AnalysisTree::Chain* treeIn = new AnalysisTree::Chain((std::vector<std::string>){fileListName}, (std::vector<std::string>){treeName});

  AnalysisTree::GenericDetector* genParticles{nullptr};

  const int ptGenFieldId = treeIn->GetConfiguration()->GetBranchConfig("Generated").GetFieldId("fGen_Pt");
  const int promptnessFieldId = treeIn->GetConfiguration()->GetBranchConfig("Generated").GetFieldId("fGen_OriginMcGen");

  treeIn->SetBranchAddress("Generated.", &genParticles);

  // TH1D* histoPtGen = new TH1D("histoPtGen", "", 2000, 0, 20);

  // ================ custom binning ================================
  std::vector<double> binEdges;
  double pt{0.};
  while(pt <= 20) {
    binEdges.emplace_back(pt);
    if(pt < 10) pt += 0.01;
    else        pt += 0.1;
  }
  TH1D* histoPtGenPrompt = new TH1D("histoPtGenPrompt", "histoPtGenPrompt", binEdges.size()-1, binEdges.data());
  TH1D* histoPtGenNonPrompt = new TH1D("histoPtGenNonPrompt", "histoPtGenNonPrompt", binEdges.size()-1, binEdges.data());
  // ================================================================

  for(const auto& h : {histoPtGenPrompt, histoPtGenNonPrompt}) {
    h->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    h->GetYaxis()->SetTitle("Entries");
  }

  const int nEntries = treeIn->GetEntries();
  for(int iEntry=0; iEntry<nEntries; ++iEntry) {
    treeIn->GetEntry(iEntry);
    for(const auto& genParticle : *(genParticles->GetChannels()) ) {
      const float pT = genParticle.GetField<float>(ptGenFieldId);
      const int promptness = genParticle.GetField<int>(promptnessFieldId);

      if(promptness == 1) histoPtGenPrompt->Fill(pT);
      else if(promptness == 2) histoPtGenNonPrompt->Fill(pT);
    }
  }

  histoPtGenPrompt->SaveAs("ptGenPrompt.root");
  histoPtGenNonPrompt->SaveAs("ptGenNonPrompt.root");
}
