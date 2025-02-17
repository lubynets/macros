#include "Configuration.hpp"
#include "Detector.hpp"
#include "Matching.hpp"
#include "Particle.hpp"
#include "PlainTreeFiller.hpp"
#include "TaskManager.hpp"

#include "TBranch.h"
#include "TFile.h"
#include "TTree.h"

#include <iostream>
#include <fstream>
#include <string>

struct IndexMap {
  std::string name_;
  std::string field_type_;
  short index_;
};

struct FicCarrier {
  float float_{-999.f};
  int int_{-999};
  char char_{static_cast<char>(-999)};
  short short_{static_cast<short>(-999)};
};

void SetAddressFIC(TBranch* branch, const IndexMap& imap, FicCarrier& ficc);
void SetFieldsFIC(const std::vector<IndexMap>& imap, AnalysisTree::Particle& particle, const std::vector<FicCarrier>& ficc);
std::vector<std::string> GetDFNames(const std::string& fileName);
bool string_to_bool(const std::string& str);

void AliceTree2AT(const std::string& fileName, bool isMC, bool isDoPlain, int maxEntries) {

  std::vector<std::string> fields_to_ignore_{"Lite_fChi2PCA",
                                             "Lite_fCpa",
                                             "Lite_fCpaXY",
                                             "Lite_fCt",
                                             "Lite_fDecayLength",
                                             "Lite_fDecayLengthXY",
                                             "Lite_fEta",
                                             "Lite_fImpactParameter0",
                                             "Lite_fImpactParameter1",
                                             "Lite_fImpactParameter2",
                                             "Lite_fM",
                                             "Lite_fMassKPi",
                                             "Lite_fPhi",
                                             "Lite_fPt",
                                             "Lite_fPtProng0",
                                             "Lite_fPtProng1",
                                             "Lite_fPtProng2"};

  bool isConfigInitialized{false};
  AnalysisTree::Configuration config_;

  AnalysisTree::Particles* candidates_{nullptr};
  AnalysisTree::BranchConfig CandidatesConfig("Candidates", AnalysisTree::DetType::kParticle);
  std::vector<IndexMap> candidateMap;
  int kfLiteSepar;
  std::vector<FicCarrier> candValues;

  AnalysisTree::Particles* simulated_{nullptr}; // MC matched with reco
  AnalysisTree::BranchConfig SimulatedConfig("Simulated", AnalysisTree::DetType::kParticle);
  std::vector<IndexMap> simulatedMap;
  std::vector<FicCarrier> simValues;

  AnalysisTree::Particles* generated_{nullptr}; // MC all decaying by 3-prong channel
  AnalysisTree::BranchConfig GeneratedConfig("Generated", AnalysisTree::DetType::kParticle);
  std::vector<IndexMap> generatedMap;
  std::vector<FicCarrier> genValues;

  AnalysisTree::Matching* cand2sim_{nullptr};

  auto CreateConfiguration = [&] (TTree* t,
                                 const std::string& prefix,
                                 AnalysisTree::BranchConfig& branch_config,
                                 std::vector<IndexMap>& vmap
                             ) {
    auto lol = t->GetListOfLeaves();
    const int nLeaves = lol->GetEntries();
    for(int iLeave=0; iLeave<nLeaves; iLeave++) {
      auto leave = lol->At(iLeave);
      const std::string fieldName = leave->GetName();
      const std::string fieldType = leave->ClassName();
      if (std::find(fields_to_ignore_.begin(), fields_to_ignore_.end(), prefix + fieldName) != fields_to_ignore_.end()) continue;
      if (fieldType == "TLeafF") {
        branch_config.AddField<float>((prefix + fieldName).c_str());
      } else if (fieldType == "TLeafI" || fieldType == "TLeafB" || fieldType == "TLeafS") {
        branch_config.AddField<int>((prefix + fieldName).c_str());
      }
      vmap.emplace_back((IndexMap){fieldName, fieldType, branch_config.GetFieldId((prefix + fieldName).c_str())});
    }
  };

  int sb_status_field_id;
  int iGlobalEntry{0};

  auto dirNames = GetDFNames(fileName);
  std::cout << "dirNames.size() = " << dirNames.size() << "\n";

  std::string fileOutName = "AnalysisTree.root";
  TFile* out_file_ = new TFile(fileOutName.c_str(), "recreate");
  TTree* tree_ = new TTree("aTree", "Analysis Tree");
  tree_->SetAutoSave(0);

  for(auto& dirname : dirNames) {
    TFile* fileIn = TFile::Open(fileName.c_str(), "read");

    TTree* treeKF = fileIn->Get<TTree>((dirname + "/O2hfcandlckf").c_str());
    TTree* treeLite = fileIn->Get<TTree>((dirname + "/O2hfcandlclite").c_str());
    TTree* treeMC = isMC ? fileIn->Get<TTree>((dirname + "/O2hfcandlcmc").c_str()) : nullptr;
    TTree* treeGen = isMC ? fileIn->Get<TTree>((dirname + "/O2hfcandlcfullp").c_str()) : nullptr;

    if(!isConfigInitialized) {
      CreateConfiguration(treeKF, "KF_", CandidatesConfig, candidateMap);
      kfLiteSepar = candidateMap.size();
      CreateConfiguration(treeLite, "Lite_", CandidatesConfig, candidateMap);
      candValues.resize(candidateMap.size());
      sb_status_field_id = std::distance(candidateMap.begin(),
                                         std::find_if(candidateMap.begin(), candidateMap.end(),
                                                      [](const IndexMap& p) { return p.name_ == "fSigBgStatus"; }));
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

        CreateConfiguration(treeGen, "Gen_", GeneratedConfig, generatedMap);
        genValues.resize(generatedMap.size());
        config_.AddBranchConfig(GeneratedConfig);
        generated_ = new AnalysisTree::Particles(GeneratedConfig.GetId());
        tree_->Branch((GeneratedConfig.GetName() + ".").c_str(), "AnalysisTree::Particles", &generated_);
      }
      config_.Print();
      isConfigInitialized = true;
    }

    candidates_->ClearChannels();
    if(isMC) {
      simulated_->ClearChannels();
      cand2sim_->Clear();

      generated_->ClearChannels();
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
      for(int iV=0; iV<genValues.size(); iV++) {
        TBranch* branch = treeGen->GetBranch(generatedMap.at(iV).name_.c_str());
        SetAddressFIC(branch, generatedMap.at(iV), genValues.at(iV));
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

      if(isMC && (candValues.at(sb_status_field_id).int_ == 1 || candValues.at(sb_status_field_id).int_ == 2)) {

        auto& simulated = simulated_->AddChannel(config_.GetBranchConfig(simulated_->GetId()));
        SetFieldsFIC(simulatedMap, simulated, simValues);

        cand2sim_->AddMatch(candidate.GetId(), simulated.GetId());
      }
      ++iGlobalEntry;
    }
    if(isMC) {
      const int nGenEntries = treeGen->GetEntries();
      for(int iEntry=0; iEntry<nGenEntries; iEntry++) {
        treeGen->GetEntry(iEntry);
        auto& generated = generated_->AddChannel(config_.GetBranchConfig(generated_->GetId()));
        SetFieldsFIC(generatedMap, generated, genValues);
      }
    }
    tree_->Fill();
    fileIn->Close();
  }

  out_file_->cd();
  config_.Write("Configuration");
  tree_->Write();
  out_file_->Close();

  if (isDoPlain) {
    std::ofstream filelist;
    filelist.open("filelist.txt");
    filelist << fileOutName + "\n";
    filelist.close();

    auto* tree_task = new AnalysisTree::PlainTreeFiller();
    tree_task->SetOutputName("PlainTree.root", "pTree");
    std::string branchname_rec = "Candidates";
    tree_task->SetInputBranchNames({branchname_rec});
    tree_task->AddBranch(branchname_rec);
    tree_task->SetIsIgnoreDefaultFields();

    auto* man = AnalysisTree::TaskManager::GetInstance();
    man->AddTask(tree_task);

    man->Init({"filelist.txt"}, {"aTree"});
    man->Run(-1);// -1 = all events
    man->Finish();
  } // isDoPlain
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./alicetree2at fileName (isMC=true isDoPlain=false nEntries=ALL)" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string fileName = argv[1];
  const bool isMC = argc > 2 ? string_to_bool(argv[2]) : true;
  const bool isDoPlain = argc > 3 ? string_to_bool(argv[3]) : false;
  const int nEntries = argc > 4 ? atoi(argv[4]) : -1;
  AliceTree2AT(fileName, isMC, isDoPlain, nEntries);

  return 0;
}

void SetAddressFIC(TBranch* branch, const IndexMap& imap, FicCarrier& ficc) {
  if     (imap.field_type_ == "TLeafF") branch->SetAddress(&ficc.float_);
  else if(imap.field_type_ == "TLeafI") branch->SetAddress(&ficc.int_);
  else if(imap.field_type_ == "TLeafB") branch->SetAddress(&ficc.char_);
  else if(imap.field_type_ == "TLeafS") branch->SetAddress(&ficc.short_);
}

void SetFieldsFIC(const std::vector<IndexMap>& imap, AnalysisTree::Particle& particle, const std::vector<FicCarrier>& ficc) {
  for(int iV=0; iV<ficc.size(); iV++) {
    if     (imap.at(iV).field_type_ == "TLeafF") particle.SetField(ficc.at(iV).float_, imap.at(iV).index_);
    else if(imap.at(iV).field_type_ == "TLeafI") particle.SetField(ficc.at(iV).int_, imap.at(iV).index_);
    else if(imap.at(iV).field_type_ == "TLeafB") particle.SetField(static_cast<int>(ficc.at(iV).char_), imap.at(iV).index_);
    else if(imap.at(iV).field_type_ == "TLeafS") particle.SetField(static_cast<int>(ficc.at(iV).short_), imap.at(iV).index_);
  }
}

std::vector<std::string> GetDFNames(const std::string& fileName) {
  TFile* fileIn = TFile::Open(fileName.c_str(), "read");
  if(fileIn == nullptr) {
    throw std::runtime_error("fileIn == nullptr");
  }

  std::vector<std::string> result;
  auto lok = fileIn->GetListOfKeys();
  for(const auto& k : *lok) {
    const std::string dirname = k->GetName();
    if(dirname == "parentFiles") continue;
    result.emplace_back(dirname);
  }
  fileIn->Close();

  return result;
}

bool string_to_bool(const std::string& str) {
  if(str == "true") return true;
  else if(str == "false") return false;
  else throw std::runtime_error("string_to_bool(): argument must be either true or false");
}
