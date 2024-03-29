void run_effmap_build(const TString simfilename, const TString recofilename)
{
  const TString y_name = "y_{CM}";
  const TString pT_name = "pT";
  const TString phi_name = "#varphi, rad";
  const int y_nbins = 22;
  const int pT_nbins = 25;
  const int phi_nbins = 30;
  const float y_beam = 1.62179;
  const float y_low = y_beam - 0.8;
  const float y_up = y_beam + 1.4;
  const float pT_low = 0.;
  const float pT_up = 2.5;
  const float phi_low = -TMath::Pi();
  const float phi_up = TMath::Pi();
  
  TFile* sim_file = TFile::Open(simfilename);
  TTree* sim_tree = sim_file->Get<TTree>("aTree");
  TFile* reco_file = TFile::Open(recofilename);
  TTree* reco_tree = reco_file->Get<TTree>("sTree");
  
  auto* sim_tracks = new AnalysisTree::Particles();
  auto* reco_tracks = new AnalysisTree::Particles();
  
  sim_tree -> SetBranchAddress("SimParticles", &sim_tracks);
  reco_tree -> SetBranchAddress("LambdaCandidates", &reco_tracks);
  
  TH2F hsim_y_pt("hsim_y_pt", "", y_nbins, y_low, y_up, pT_nbins, pT_low, pT_up);
  TH2F hrec_y_pt("hrec_y_pt", "", y_nbins, y_low, y_up, pT_nbins, pT_low, pT_up);
  TH2F hsim_phi_pt("hsim_phi_pt", "", phi_nbins, phi_low, phi_up, pT_nbins, pT_low, pT_up);
  TH2F hrec_phi_pt("hrec_phi_pt", "", phi_nbins, phi_low, phi_up, pT_nbins, pT_low, pT_up);
  TH2F hsim_phi_y("hsim_phi_y", "", phi_nbins, phi_low, phi_up, y_nbins, y_low, y_up);
  TH2F hrec_phi_y("hrec_phi_y", "", phi_nbins, phi_low, phi_up, y_nbins, y_low, y_up);
  
  hsim_y_pt.GetXaxis()->SetTitle(y_name);
  hrec_y_pt.GetXaxis()->SetTitle(y_name);
  hsim_phi_pt.GetXaxis()->SetTitle(phi_name);
  hrec_phi_pt.GetXaxis()->SetTitle(phi_name);
  hsim_phi_y.GetXaxis()->SetTitle(phi_name);
  hrec_phi_y.GetXaxis()->SetTitle(phi_name);
  
  hsim_y_pt.GetYaxis()->SetTitle(pT_name);
  hrec_y_pt.GetYaxis()->SetTitle(pT_name);
  hsim_phi_pt.GetYaxis()->SetTitle(pT_name);
  hrec_phi_pt.GetYaxis()->SetTitle(pT_name);
  hsim_phi_y.GetYaxis()->SetTitle(y_name);
  hrec_phi_y.GetYaxis()->SetTitle(y_name);
  
  const int sim_entries = sim_tree->GetEntries();
  for(int iEvent=0; iEvent<sim_entries; iEvent++)
  {
    sim_tree->GetEntry(iEvent);
    for(const auto& simtrack : *(sim_tracks->GetChannels()))
    {
      if(simtrack.GetPid() != 3122) continue;
      
      hsim_y_pt.Fill(simtrack.GetRapidity(), simtrack.GetPt());
      hsim_phi_pt.Fill(simtrack.GetPhi(), simtrack.GetPt());
      hsim_phi_y.Fill(simtrack.GetPhi(), simtrack.GetRapidity());
    }    
  }
  
  const int rec_entries = reco_tree->GetEntries();
  for(int iEvent=0; iEvent<rec_entries; iEvent++)
  {
    reco_tree->GetEntry(iEvent);
    for(const auto& recotrack : *(reco_tracks->GetChannels()))
    {
      if(recotrack.GetField<int>(5) == 0) continue;
      
      hrec_y_pt.Fill(recotrack.GetRapidity(), recotrack.GetPt());
      hrec_phi_pt.Fill(recotrack.GetPhi(), recotrack.GetPt());
      hrec_phi_y.Fill(recotrack.GetPhi(), recotrack.GetRapidity());
    }
  }
  
  TFile fileOut("map_pt_y_phi.root", "recreate");
  fileOut.cd();
  
  hsim_y_pt.Write();
  hrec_y_pt.Write();
  hsim_phi_pt.Write();
  hrec_phi_pt.Write();
  hsim_phi_y.Write();
  hrec_phi_y.Write();
  
  fileOut.Close();  
}
  