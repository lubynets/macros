void countNEvents() {
  const std::string filePath = "/lustre/alice/users/lubynets/ali2atree/outputs/data/lhc22.apass7/noConstr/set1/";
  const int fileFrom = 1;
  const int fileTo = 976;


  long long int nEntries{0};
  for(int iFile=fileFrom; iFile<=fileTo; iFile++) {
    TFile* fileIn = TFile::Open((filePath + "AnalysisTree." + std::to_string(iFile) + ".root").c_str(), "read");
    TTree* treeIn = fileIn->Get<TTree>("aTree");
    nEntries += treeIn->GetEntries();
    fileIn->Close();
  }
  std::cout << "In files from " << std::to_string(fileFrom) << " to " << std::to_string(fileTo) << " nEntries = " << nEntries << "\n";
}
