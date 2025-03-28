void ListVarsToSkipWriter(const std::string& fileName, bool eachNewLine=false) {
  const std::vector<std::string> varsToPreserve {
    "Candidates_KF_fChi2PrimKaon",
    "Candidates_KF_fChi2PrimPion",
    "Candidates_KF_fChi2PrimProton",
    "Candidates_KF_fChi2GeoPionKaon",
    "Candidates_KF_fChi2GeoProtonKaon",
    "Candidates_KF_fChi2GeoProtonPion",
    "Candidates_KF_fChi2Geo",
    "Candidates_KF_fDecayLengthNormalised",
    "Candidates_KF_fChi2Topo",
    "Candidates_KF_fSigBgStatus",
    "Candidates_KF_fMassInv",
    "Candidates_Lite_fPt",
    "Candidates_Lite_fT"
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
