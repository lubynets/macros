void AO2DFileReducer() {
  std::string fileInName = "/lustre/alice/users/lubynets/ao2ds/data/2022/LHC22o/526641/apass7/0630/o2_ctf_run00526641_orbit0206830848_tf0000000001_epn160/001/AO2D.root";
  std::vector<std::string> dataFrameNames{"DF_2261906078563584", "DF_2261906078563840"};
  std::string fileOutName = "fileOut.root";

  TFile* fileIn = TFile::Open(fileInName.c_str(), "open");
  TFile* fileOut = TFile::Open(fileOutName.c_str(), "recreate");

  TMap* tmap = (TMap*)fileIn->Get<TMap>("metaData")->Clone();
  fileOut->cd();
  tmap->Write("metaData", 1);

  for(auto& dFN : dataFrameNames) {
    TDirectoryFile* tdf = fileIn->Get<TDirectoryFile>(dFN.c_str());
    TDirectoryFile* newtdf = new TDirectoryFile(dFN.c_str(), dFN.c_str());
    fileOut->cd();
    newtdf->Write();
    newtdf->cd();
    TList* lok = tdf->GetListOfKeys();
    for(const auto& k : *lok) {
      std::string treeName = k->GetName();
      TTree* tree = (TTree*)fileIn->Get<TTree>((dFN + "/" + treeName).c_str())->CloneTree();
      tree->Write();
    }
    fileOut->cd();
  }

  fileOut->Close();
  fileIn->Close();
}
