struct IndexMap {
  std::string name_;
  std::string field_type_;
  short index_;
};

struct FicCarrier {
  float float_;
  int int_;
  char char_;
};

void SetAddressFIC(TBranch* branch, const IndexMap& imap, FicCarrier& ficc);
void SetFieldsFIC(const std::vector<IndexMap>& imap, AnalysisTree::Particle& particle, const std::vector<FicCarrier>& ficc);

void AliceTree2AT(const std::string& fileName, bool isMC=true, int maxEntries=-1) {
  TFile* fileIn = TFile::Open(fileName.c_str());
  if(fileIn == nullptr) {
    throw std::runtime_error("fileIn == nullptr");
  }

  bool isConfigInitialized{false};
  AnalysisTree::Configuration config_;

  AnalysisTree::Particles* candidates_{nullptr};
  AnalysisTree::BranchConfig CandidatesConfig("Candidates", AnalysisTree::DetType::kParticle);
  std::vector<IndexMap> candidateMap;
  int kfLiteSepar;
  std::vector<FicCarrier> candValues;

  AnalysisTree::Particles* simulated_{nullptr};
  AnalysisTree::BranchConfig SimulatedConfig("Simulated", AnalysisTree::DetType::kParticle);
  std::vector<IndexMap> simulatedMap;
  std::vector<FicCarrier> simValues;

  AnalysisTree::Matching* cand2sim_{nullptr};

  TFile* out_file_ = new TFile("AnalysisTree.root", "recreate");

  TTree* tree_{nullptr};

  auto CreateConfiguration = [&] (TTree* t,
                                  std::string prefix,
                                  AnalysisTree::BranchConfig& branch_config,
                                  std::vector<IndexMap>& vmap
                                  ) {
    auto lol = t->GetListOfLeaves();
    const int nLeaves = lol->GetEntries();
    for(int iLeave=0; iLeave<nLeaves; iLeave++) {
      auto leave = lol->At(iLeave);
      const std::string fieldName = leave->GetName();
      const std::string fieldType = leave->ClassName();
      if (fieldType == "TLeafF") {
        branch_config.AddField<float>((prefix + fieldName).c_str());
      } else if (fieldType == "TLeafI" || fieldType == "TLeafB") {
        branch_config.AddField<int>((prefix + fieldName).c_str());
      }
      vmap.emplace_back((IndexMap){fieldName, fieldType, branch_config.GetFieldId((prefix + fieldName).c_str())});
    }
  };

  int sb_status_field_id;
  int iGlobalEntry{0};

  auto lok = fileIn->GetListOfKeys();
  for(const auto& k : *lok) {
    const std::string dirname = k->GetName();
    if(dirname == "parentFiles") continue;

    TTree* treeKF = fileIn->Get<TTree>((dirname + "/O2hfcandlckf").c_str());
    TTree* treeLite = fileIn->Get<TTree>((dirname + "/O2hfcandlclite").c_str());
    TTree* treeMC = isMC ? fileIn->Get<TTree>((dirname + "/O2hfcandlcmc").c_str()) : nullptr;

    if(!isConfigInitialized) {
      tree_ = new TTree("aTree", "Analysis Tree");

      CreateConfiguration(treeKF, "KF_", CandidatesConfig, candidateMap);
      kfLiteSepar = candidateMap.size();
      CreateConfiguration(treeLite, "Lite_", CandidatesConfig, candidateMap);
      candValues.resize(candidateMap.size());
      sb_status_field_id = std::find_if(candidateMap.begin(), candidateMap.end(),
                                        [](const IndexMap& p) { return p.name_ == "fSigBgStatus"; })->index_;
      config_.AddBranchConfig(CandidatesConfig);
      candidates_ = new AnalysisTree::Particles(CandidatesConfig.GetId());
      tree_->Branch((CandidatesConfig.GetName() + ".").c_str(), "AnalysisTree::Particles", &candidates_);
      if(isMC) {
        CreateConfiguration(treeMC, "Sim_", SimulatedConfig, simulatedMap);
        simValues.resize(simulatedMap.size());
        config_.AddBranchConfig(SimulatedConfig);
        simulated_ = new AnalysisTree::Particles(SimulatedConfig.GetId());
        tree_->Branch((SimulatedConfig.GetName() + ".").c_str(), "AnalysisTree::Particles", &simulated_);
        cand2sim_ = new AnalysisTree::Matching(CandidatesConfig.GetId(), SimulatedConfig.GetId());
        config_.AddMatch(cand2sim_);
        tree_->Branch((CandidatesConfig.GetName() + "2" + SimulatedConfig.GetName() + ".").c_str(), "AnalysisTree::Matching", &cand2sim_);
      }
      isConfigInitialized = true;
    }

    candidates_->ClearChannels();
    if(isMC) {
      simulated_->ClearChannels();
      cand2sim_->Clear();
    }

    for(int iV=0; iV<candValues.size(); iV++) {
      auto treeRec = iV<kfLiteSepar ? treeKF : treeLite;
      TBranch* branch = treeRec->GetBranch(candidateMap.at(iV).name_.c_str());
      SetAddressFIC(branch, candidateMap.at(iV), candValues.at(iV));
    }
    if(isMC) {
      for(int iV=0; iV<simValues.size(); iV++) {
        TBranch* branch = treeMC->GetBranch(simulatedMap.at(iV).name_.c_str());
        SetAddressFIC(branch, simulatedMap.at(iV), simValues.at(iV));
      }
    }

    const int nEntries = treeKF->GetEntries();
    if(treeLite->GetEntries() != nEntries || (isMC && treeMC->GetEntries() != nEntries)) {
      throw std::runtime_error("treeLite->GetEntries() != nEntries || treeMC->GetEntries() != nEntries");
    }

    for(int iEntry=0; iEntry<nEntries; iEntry++) {
      if(maxEntries > 0 && iGlobalEntry >= maxEntries) break;
      treeKF->GetEntry(iEntry);
      treeLite->GetEntry(iEntry);
      if(isMC) treeMC->GetEntry(iEntry);

      auto& candidate = candidates_->AddChannel(config_.GetBranchConfig(candidates_->GetId()));
      SetFieldsFIC(candidateMap, candidate, candValues);

      if(isMC && (candValues.at(sb_status_field_id).int_ == 0 || candValues.at(sb_status_field_id).int_ == 1)) {

        auto& simulated = simulated_->AddChannel(config_.GetBranchConfig(simulated_->GetId()));
        SetFieldsFIC(simulatedMap, simulated, simValues);

        cand2sim_->AddMatch(candidate.GetId(), simulated.GetId());
      }
      ++iGlobalEntry;
    }
    tree_->Fill();
  }

  out_file_->cd();
  config_.Write("Configuration");
  tree_->Write();
  out_file_->Close();

  fileIn->Close();
}

void SetAddressFIC(TBranch* branch, const IndexMap& imap, FicCarrier& ficc) {
  if     (imap.field_type_ == "TLeafF") branch->SetAddress(&ficc.float_);
  else if(imap.field_type_ == "TLeafI") branch->SetAddress(&ficc.int_);
  else if(imap.field_type_ == "TLeafB") branch->SetAddress(&ficc.char_);
}

void SetFieldsFIC(const std::vector<IndexMap>& imap, AnalysisTree::Particle& particle, const std::vector<FicCarrier>& ficc) {
  for(int iV=0; iV<ficc.size(); iV++) {
    if     (imap.at(iV).field_type_ == "TLeafF") particle.SetField(ficc.at(iV).float_, imap.at(iV).index_);
    else if(imap.at(iV).field_type_ == "TLeafI") particle.SetField(ficc.at(iV).int_, imap.at(iV).index_);
    else if(imap.at(iV).field_type_ == "TLeafB") particle.SetField(static_cast<int>(ficc.at(iV).char_), imap.at(iV).index_);
  }
}
