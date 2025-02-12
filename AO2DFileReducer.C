void AO2DFileReducer() {
  std::string fileInName = "/lustre/alice/users/lubynets/CSTlc/outputs/mc/lhc24e3_tm/all/noConstr/AnalysisResults_trees.1.root";
  std::vector<std::string> dataFrameNames{"DF_2261906081273110", "DF_2261906080459030", "DF_2261906080459158", "DF_2261906080459286", "DF_2261906080459414"};
  std::string fileOutName = "fileOut.root";

  TFile* fileIn = TFile::Open(fileInName.c_str(), "open");
  TFile* fileOut = TFile::Open(fileOutName.c_str(), "recreate");

//   TMap* tmap = (TMap*)fileIn->Get<TMap>("metaData")->Clone();
//   fileOut->cd();
//   tmap->Write("metaData", 1);

  for(auto& dFN : dataFrameNames) {
    TDirectoryFile* tdf = fileIn->Get<TDirectoryFile>(dFN.c_str());
    TDirectoryFile* newtdf = new TDirectoryFile(dFN.c_str(), dFN.c_str());
    fileOut->cd();
    newtdf->Write();
    newtdf->cd();
    TList* lok = tdf->GetListOfKeys();
    for(const auto& k : *lok) {
      std::string treeName = k->GetName();
      std::cout << "treeName = " << treeName << "\n";
      TTree* tree = (TTree*)fileIn->Get<TTree>((dFN + "/" + treeName).c_str())->CloneTree();
      tree->Write();
    }
    fileOut->cd();
  }

  fileOut->Close();
  fileIn->Close();
}
