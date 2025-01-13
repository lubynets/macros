void AliceTree2AT(const std::string& fileName, bool isMC=true, int nEntries=-1) {
  TFile* fileIn = TFile::Open(fileName.c_str());
  if(fileIn == nullptr) {
    throw std::runtime_error("fileIn == nullptr");
  }

  bool isConfigInitialized{false};
  AnalysisTree::Configuration config_;

  AnalysisTree::Particles* candidates_{nullptr};
  AnalysisTree::BranchConfig CandidatesConfig("Candidates", AnalysisTree::DetType::kParticle);
  std::vector<std::pair<int, std::string>> candidateMapFloats;
  std::vector<std::pair<int, std::string>> candidateMapInts;
  int kfLiteSeparF;
  int kfLiteSeparI;
  std::vector<float> candValuesF;
  std::vector<int> candValuesI;

  AnalysisTree::Particles* simulated_{nullptr};
  AnalysisTree::BranchConfig SimulatedConfig("Simulated", AnalysisTree::DetType::kParticle);
  std::vector<std::pair<int, std::string>> simulatedMapFloats;
  std::vector<std::pair<int, std::string>> simulatedMapInts;
  std::vector<float> simValuesF;
  std::vector<int> simValuesI;

  AnalysisTree::Matching* cand2sim_{nullptr};

  TFile* out_file_ = new TFile("AnalysisTree.root", "recreate");

  TTree* tree_{nullptr};

  auto CreateConfiguration = [&] (TTree* t,
                                  std::string prefix,
                                  AnalysisTree::BranchConfig& branch_config,
                                  std::vector<std::pair<int, std::string>>& mapF,
                                  std::vector<std::pair<int, std::string>>& mapI) {
    auto lol = t->GetListOfLeaves();
    const int nLeaves = lol->GetEntries();
    for(int iLeave=0; iLeave<nLeaves; iLeave++) {
      auto leave = lol->At(iLeave);
      const std::string fieldName = leave->GetName();
      const std::string fieldType = leave->ClassName();
      if (fieldType == "TLeafF") {
        branch_config.AddField<float>((prefix + fieldName).c_str());
        mapF.emplace_back(std::make_pair(branch_config.GetFieldId((prefix + fieldName).c_str()), fieldName));
      } else if (fieldType == "TLeafI") {
        branch_config.AddField<int>((prefix + fieldName).c_str());
        mapI.emplace_back(std::make_pair(branch_config.GetFieldId((prefix + fieldName).c_str()), fieldName));
      }
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

      CreateConfiguration(treeKF, "KF_", CandidatesConfig, candidateMapFloats, candidateMapInts);
      kfLiteSeparF = candidateMapFloats.size();
      kfLiteSeparI = candidateMapInts.size();
      CreateConfiguration(treeLite, "Lite_", CandidatesConfig, candidateMapFloats, candidateMapInts);
      candValuesF.resize(candidateMapFloats.size());
      candValuesI.resize(candidateMapInts.size());
      sb_status_field_id = std::find_if(candidateMapInts.begin(), candidateMapInts.end(),
                                        [](const std::pair<int, std::string>& p) { return p.second == "fSigBgStatus"; })->first;
      config_.AddBranchConfig(CandidatesConfig);
      candidates_ = new AnalysisTree::Particles(CandidatesConfig.GetId());
      tree_->Branch((CandidatesConfig.GetName() + ".").c_str(), "AnalysisTree::Particles", &candidates_);
      if(isMC) {
        CreateConfiguration(treeMC, "Sim_", SimulatedConfig, simulatedMapFloats, simulatedMapInts);
        simValuesF.resize(simulatedMapFloats.size());
        simValuesI.resize(simulatedMapInts.size());
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

    for(int iV=0; iV<candValuesF.size(); iV++) {
      if(iV<kfLiteSeparF) treeKF->SetBranchAddress(candidateMapFloats.at(iV).second.c_str(), &candValuesF.at(iV));
      else                treeLite->SetBranchAddress(candidateMapFloats.at(iV).second.c_str(), &candValuesF.at(iV));
    }
    for(int iV=0; iV<candValuesI.size(); iV++) {
      if(iV<kfLiteSeparI) treeKF->SetBranchAddress(candidateMapInts.at(iV).second.c_str(), &candValuesI.at(iV));
      else                treeLite->SetBranchAddress(candidateMapInts.at(iV).second.c_str(), &candValuesI.at(iV));
    }
    if(isMC) {
      for(int iV=0; iV<std::max(candValuesF.size(), candValuesI.size()); iV++) {
        if(iV<simValuesF.size()) treeMC->SetBranchAddress(simulatedMapFloats.at(iV).second.c_str(), &simValuesF.at(iV));
        if(iV<simValuesI.size()) treeMC->SetBranchAddress(simulatedMapInts.at(iV).second.c_str(), &simValuesI.at(iV));
      }
    }

    const int nEntries = treeKF->GetEntries();
    if(treeLite->GetEntries() != nEntries || (isMC && treeMC->GetEntries() != nEntries)) {
      throw std::runtime_error("treeLite->GetEntries() != nEntries || treeMC->GetEntries() != nEntries");
    }

    for(int iEntry=0; iEntry<nEntries; iEntry++) {
      if(nEntries > 0 && iGlobalEntry >= nEntries) break;
      treeKF->GetEntry(iEntry);
      treeLite->GetEntry(iEntry);
      if(isMC) treeMC->GetEntry(iEntry);

      auto& candidate = candidates_->AddChannel(config_.GetBranchConfig(candidates_->GetId()));

      for(int iV=0; iV<std::max(candValuesF.size(), candValuesI.size()); iV++) {
        if(iV<candValuesF.size()) candidate.SetField(candValuesF.at(iV), candidateMapFloats.at(iV).first);
        if(iV<candValuesI.size()) candidate.SetField(candValuesI.at(iV), candidateMapInts.at(iV).first);
      }

      if(isMC && candValuesI.at(sb_status_field_id) > 0) {
        auto& simulated = simulated_->AddChannel(config_.GetBranchConfig(simulated_->GetId()));

        for(int iV=0; iV<std::max(simValuesF.size(), simValuesI.size()); iV++) {
          if(iV<simValuesF.size()) simulated.SetField(simValuesF.at(iV), simulatedMapFloats.at(iV).first);
          if(iV<simValuesI.size()) simulated.SetField(simValuesI.at(iV), simulatedMapInts.at(iV).first);
        }
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
