void ListVarsToSkipWriter(const std::string& fileName, bool eachNewLine=false) {
  const std::vector<std::string> varsToPreserve {
    "fLiteImpactParameter0",
    "fLiteImpactParameter1",
    "fLiteImpactParameter2",
    "fLiteChi2PCA",
    "fLiteCpa",
    "fLiteCpaXY",
    "fLiteDecayLengthXY",
    "fLiteDecayLength",
    "fKFPt",
    "fKFMassInv",
    "fKFT",
    "fLiteY",
    "fKFSigBgStatus"
  };

  TFile* fileIn = TFile::Open(fileName.c_str(), "open");
  TTree* treeIn = fileIn->Get<TTree>("pTree");

  std::vector<std::string> varsToSkip;

  for(int iLeaf=0, nLeaves=treeIn->GetListOfLeaves()->GetEntries(); iLeaf<nLeaves; iLeaf++) {
    const std::string varName = treeIn->GetListOfLeaves()->At(iLeaf)->GetName();
    if(std::find(varsToPreserve.begin(), varsToPreserve.end(), varName) == varsToPreserve.end()) {
      varsToSkip.emplace_back(varName);
    }
  }
  fileIn->Close();

  const std::string separ = eachNewLine ? "\n" : " ";

  for(auto& vTS : varsToSkip) {
    std::cout << "\'" << vTS << "\'," + separ;
  }
}
