#include "Configuration.hpp"
#include "Detector.hpp"
#include "EventHeader.hpp"
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
template <typename T>
void SetFieldsFIC(const std::vector<IndexMap>& imap, T& obj, const std::vector<FicCarrier>& ficc);
std::vector<std::string> GetDFNames(const std::string& fileName);
int DetermineFieldIdByName(const std::vector<IndexMap>& iMap, const std::string& name);
bool string_to_bool(const std::string& str);
std::vector<int> findPositions(const std::vector<int>& vec, int M);

void AliceTree2AT(const std::string& fileName, bool isMC, bool isDoPlain, int maxEntries) {

  const std::vector<std::string> fields_to_ignore_{
    "fLiteChi2PCA",
    "fLiteCpa",
    "fLiteCpaXY",
    "fLiteCt",
    "fLiteDecayLength",
    "fLiteDecayLengthXY",
    "fLiteEta",
    "fLiteImpactParameter0",
    "fLiteImpactParameter1",
    "fLiteImpactParameter2",
    "fLiteM",
    "fLiteMassKPi",
    "fLitePhi",
    "fLitePt",
    "fLitePtProng0",
    "fLitePtProng1",
    "fLitePtProng2"
  };

  const std::vector<std::string> fields_to_preserve_ {};

  if(!fields_to_ignore_.empty() && !fields_to_preserve_.empty()) throw std::runtime_error("!fields_to_ignore_.empty() && !fields_to_preserve_.empty()");

  bool isConfigInitialized{false};
  AnalysisTree::Configuration config_;

  AnalysisTree::EventHeader* eve_header_{nullptr};
  AnalysisTree::BranchConfig EventsConfig("Events", AnalysisTree::DetType::kEventHeader);
  std::vector<IndexMap> eventsMap;
  std::vector<FicCarrier> eventValues;

  AnalysisTree::Particles* candidates_{nullptr};
  AnalysisTree::BranchConfig CandidatesConfig("Candidates", AnalysisTree::DetType::kParticle);
  std::vector<IndexMap> candidateMap;
  int kfLiteSepar, liteCollIdSepar;
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
      const std::string prefixedFieldName = "f" + prefix + fieldName.substr(1, fieldName.size());
      const std::string fieldType = leave->ClassName();
      if (!fields_to_ignore_.empty() && (std::find(fields_to_ignore_.begin(), fields_to_ignore_.end(), prefixedFieldName) != fields_to_ignore_.end())) continue;
      if (!fields_to_preserve_.empty() && (std::find(fields_to_preserve_.begin(), fields_to_preserve_.end(), prefixedFieldName) == fields_to_preserve_.end())) continue;
      if (fieldType == "TLeafF") {
        branch_config.AddField<float>(prefixedFieldName);
      } else if (fieldType == "TLeafI" || fieldType == "TLeafB" || fieldType == "TLeafS") {
        branch_config.AddField<int>(prefixedFieldName);
      }
      vmap.emplace_back((IndexMap){fieldName, fieldType, branch_config.GetFieldId(prefixedFieldName)});
    }
  };

  int sb_status_field_id;
  int collision_id_field_id_in_evehead;
  int collision_id_field_id_in_cand;
  int iGlobalEntry{0};

  auto dirNames = GetDFNames(fileName);
  std::cout << "dirNames.size() = " << dirNames.size() << "\n";

  std::string fileOutName = "AnalysisTree.root";
  TFile* out_file_ = new TFile(fileOutName.c_str(), "recreate");
  TTree* tree_ = new TTree("aTree", "Analysis Tree");
  tree_->SetAutoSave(0);

  for(auto& dirname : dirNames) {
    TFile* fileIn = TFile::Open(fileName.c_str(), "read");
    bool is_gentree_processed{false};

    TTree* treeKF = fileIn->Get<TTree>((dirname + "/O2hfcandlckf").c_str());
    TTree* treeLite = fileIn->Get<TTree>((dirname + "/O2hfcandlclite").c_str());
    TTree* treeCollId = fileIn->Get<TTree>((dirname + "/O2hfcollidlclite").c_str());
    TTree* treeMC = isMC ? fileIn->Get<TTree>((dirname + "/O2hfcandlcmc").c_str()) : nullptr;
    TTree* treeGen = isMC ? fileIn->Get<TTree>((dirname + "/O2hfcandlcfullp").c_str()) : nullptr;
    TTree* treeEvent = fileIn->Get<TTree>((dirname + "/O2hfcandlcfullev").c_str());

    if(!isConfigInitialized) {
      CreateConfiguration(treeEvent, "Ev", EventsConfig, eventsMap);
      eventValues.resize(eventsMap.size());
      collision_id_field_id_in_evehead = DetermineFieldIdByName(eventsMap, "fIndexCollisions");
      config_.AddBranchConfig(EventsConfig);
      eve_header_ = new AnalysisTree::EventHeader(EventsConfig.GetId());
      eve_header_->Init(EventsConfig);
      tree_->Branch((EventsConfig.GetName() + ".").c_str(), "AnalysisTree::EventHeader", &eve_header_);

      CreateConfiguration(treeKF, "KF", CandidatesConfig, candidateMap);
      kfLiteSepar = candidateMap.size();
      CreateConfiguration(treeLite, "Lite", CandidatesConfig, candidateMap);
      liteCollIdSepar = candidateMap.size();
      CreateConfiguration(treeCollId, "", CandidatesConfig, candidateMap);
      candValues.resize(candidateMap.size());
      sb_status_field_id = DetermineFieldIdByName(candidateMap, "fSigBgStatus");
      collision_id_field_id_in_cand = DetermineFieldIdByName(candidateMap, "fIndexCollisions");
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

    for(int iV=0; iV<eventValues.size(); iV++) {
      TBranch* branch = treeEvent->GetBranch(eventsMap.at(iV).name_.c_str());
      SetAddressFIC(branch, eventsMap.at(iV), eventValues.at(iV));
    }

    for(int iV=0; iV<candValues.size(); iV++) {
      auto treeRec = iV<kfLiteSepar ? treeKF : iV<liteCollIdSepar ? treeLite : treeCollId;
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

    const int nEntriesKF = treeKF->GetEntries();
    if(treeLite->GetEntries() != nEntriesKF || treeCollId->GetEntries() != nEntriesKF || (isMC && treeMC->GetEntries() != nEntriesKF)) {
      std::cout << "treeLite->GetEntries() = " << treeLite->GetEntries() << "\n";
      std::cout << "treeKF->GetEntries() = " << treeKF->GetEntries() << "\n";
      std::cout << "treeCollId->GetEntries() = " << treeCollId->GetEntries() << "\n";
      if(isMC) std::cout << "treeMC->GetEntries() = " << treeMC->GetEntries() << "\n";
      throw std::runtime_error("treeLite->GetEntries() != nEntriesKF || treeCollId->GetEntries() != nEntriesKF || (isMC && treeMC->GetEntries() != nEntriesKF)");
    }

    std::vector<int> candidateCollisionIndices;
    for(int iEntryKF = 0; iEntryKF<nEntriesKF; iEntryKF++) {
      treeKF->GetEntry(iEntryKF);
      treeCollId->GetEntry(iEntryKF);
      candidateCollisionIndices.emplace_back(candValues.at(collision_id_field_id_in_cand).int_);
    }

    const int nEntriesEve = treeEvent->GetEntries();
    for(int iEntryEve=0; iEntryEve<nEntriesEve && (maxEntries<0 || iGlobalEntry<maxEntries); iEntryEve++, iGlobalEntry++) {
      candidates_->ClearChannels();
      if(isMC) {
        simulated_->ClearChannels();
        cand2sim_->Clear();
        generated_->ClearChannels();
      }

      treeEvent->GetEntry(iEntryEve);
      SetFieldsFIC(eventsMap, *eve_header_, eventValues);

      const int indexCollision = eventValues.at(collision_id_field_id_in_evehead).int_;
      auto candidatesOfThisCollisionIndices = findPositions(candidateCollisionIndices, indexCollision);

      for(auto& cOTI : candidatesOfThisCollisionIndices) {
        treeKF->GetEntry(cOTI);
        treeLite->GetEntry(cOTI);
        treeCollId->GetEntry(cOTI);
        if(isMC) treeMC->GetEntry(cOTI);

        auto& candidate = candidates_->AddChannel(config_.GetBranchConfig(candidates_->GetId()));
        SetFieldsFIC(candidateMap, candidate, candValues);

        if(isMC && (candValues.at(sb_status_field_id).int_ == 1 || candValues.at(sb_status_field_id).int_ == 2)) {

          auto& simulated = simulated_->AddChannel(config_.GetBranchConfig(simulated_->GetId()));
          SetFieldsFIC(simulatedMap, simulated, simValues);

          cand2sim_->AddMatch(candidate.GetId(), simulated.GetId());
        }
      } // KF entries
      if(isMC && !is_gentree_processed) {
        const int nGenEntries = treeGen->GetEntries();
        for(int iEntry=0; iEntry<nGenEntries; iEntry++) {
          treeGen->GetEntry(iEntry);
          auto& generated = generated_->AddChannel(config_.GetBranchConfig(generated_->GetId()));
          SetFieldsFIC(generatedMap, generated, genValues);
        } // Gen entries
        is_gentree_processed = true;
      } // isMC && !is_gentree_processed
      tree_->Fill();
    } // event entries
    fileIn->Close();
  } // dirNames

  out_file_->cd();
  config_.Write("Configuration");
  tree_->Write();
  out_file_->Close();

  if (isDoPlain) {
    std::ofstream filelist;
    filelist.open("filelist.txt");
    filelist << fileOutName + "\n";
    filelist.close();

    const std::array<double, 4> sidebands{2.12, 2.20, 2.38, 2.42};
    AnalysisTree::SimpleCut invMassCut = AnalysisTree::SimpleCut({"Candidates.fKFMassInv"}, [&] (const std::vector<double>& par) { return (par[0]>sidebands.at(0) && par[0]<sidebands.at(1)) || (par[0]>sidebands.at(2) && par[0]<sidebands.at(3)); });
    AnalysisTree::Cuts* sideBandCuts = new AnalysisTree::Cuts("sideBandCuts", {invMassCut});

    AnalysisTree::SimpleCut signalCut = AnalysisTree::RangeCut("Candidates.fKFSigBgStatus", 0.9, 2.1);
    AnalysisTree::Cuts* signalCuts = new AnalysisTree::Cuts("signalCuts", {signalCut});

    auto* tree_task = new AnalysisTree::PlainTreeFiller();
    tree_task->SetOutputName("PlainTree.root", "pTree");
    std::string branchname_rec = "Candidates";
    tree_task->SetInputBranchNames({branchname_rec});
    tree_task->AddBranch(branchname_rec);
//     tree_task->AddBranchCut(signalCuts);
    tree_task->AddBranchCut(sideBandCuts);
    const std::vector<std::string> fields_to_preserve {
      "fKFChi2PrimProton",
      "fKFChi2PrimKaon",
      "fKFChi2PrimPion",
      "fKFChi2Geo",
      "fKFChi2Topo",
      "fKFDecayLengthNormalised",
      "KFNSigTpcPr",
      "KFNSigTpcKa",
      "KFNSigTpcPi",
      "fKFT",
      "fKFPt",
      "fKFMassInv",
      "fKFSigBgStatus"
    };
    tree_task->SetFieldsToPreserve(fields_to_preserve);
    tree_task->SetIsPrependLeavesWithBranchName(false);

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

template <typename T>
void SetFieldsFIC(const std::vector<IndexMap>& imap, T& obj, const std::vector<FicCarrier>& ficc) {
  for(int iV=0; iV<ficc.size(); iV++) {
    if     (imap.at(iV).field_type_ == "TLeafF") obj.SetField(ficc.at(iV).float_, imap.at(iV).index_);
    else if(imap.at(iV).field_type_ == "TLeafI") obj.SetField(ficc.at(iV).int_, imap.at(iV).index_);
    else if(imap.at(iV).field_type_ == "TLeafB") obj.SetField(static_cast<int>(ficc.at(iV).char_), imap.at(iV).index_);
    else if(imap.at(iV).field_type_ == "TLeafS") obj.SetField(static_cast<int>(ficc.at(iV).short_), imap.at(iV).index_);
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

int DetermineFieldIdByName(const std::vector<IndexMap>& iMap, const std::string& name) {
  auto distance = std::distance(iMap.begin(),std::find_if(iMap.begin(), iMap.end(), [&name](const IndexMap& p) { return p.name_ == name; }));
  if(distance == iMap.size()) throw std::runtime_error("DetermineFieldIdByName(): name " + name + " is missing");
  return distance;
}

bool string_to_bool(const std::string& str) {
  if(str == "true") return true;
  else if(str == "false") return false;
  else throw std::runtime_error("string_to_bool(): argument must be either true or false");
}

std::vector<int> findPositions(const std::vector<int>& vec, int M) {
  std::vector<int> positions;
  for (int i = 0; i < vec.size(); ++i) {
    if (vec[i] == M) {
      positions.push_back(i); // Store the index
    }
  }
  return positions;
}
