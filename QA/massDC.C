void massDC(std::string filelist_eve, std::string filelist_pfs, int pdg=3122)
{
  AnalysisTree::Chain* treeIn = new AnalysisTree::Chain(std::vector<std::string>({filelist_eve, filelist_pfs}), std::vector<std::string>({"aTree", "pTree"}));
  auto* config = treeIn->GetConfiguration();
  
  auto* reco_tracks = new AnalysisTree::Particles();
  auto* event_header = new AnalysisTree::EventHeader();
   
  treeIn->SetBranchAddress("Candidates.", &reco_tracks);
  treeIn->SetBranchAddress("RecEventHeader.", &event_header);
  
  int centrality_id_ = config->GetBranchConfig("RecEventHeader").GetFieldId("centrality_tracks");
  int is_signal_id_ = config->GetBranchConfig("Candidates").GetFieldId("generation");
  
//   const float y_beam = 1.62179;
  const float y_beam = 0.985344;
  
  int mass_nbins;
  float mass_low;
  float mass_up;

  int   y_nbins;
  float y_low;
  float y_up;
  int   pT_nbins;
  float pT_low;
  float pT_up;
  std::vector<double> C_binranges;

  if(std::fabs(y_beam - 1.62179)<0.01) {
    C_binranges = {0, 5, 10, 15, 20, 30, 40, 70};
    if(pdg == 3122) {
      mass_nbins = 1200;
      mass_low = 1.045;
      mass_up = 1.185;
      y_nbins = 6;
      y_low = y_beam-0.75;
      y_up = y_beam+1.05;
      pT_nbins = 3;
      pT_low = 0.2;
      pT_up = 1.4;
    } else if(pdg == 310) {
      mass_nbins = 1200;
      mass_low = 0.277;
      mass_up = 0.717;
      y_nbins = 5;
      y_low = y_beam-0.45;
      y_up = y_beam+1.05;
      pT_nbins = 3;
      pT_low = 0;
      pT_up = 1.5;
    } else if(pdg == 3312) {
      mass_nbins = 1200;
      mass_low = 1.251;
      mass_up = 1.391;
      y_nbins = 5;
      y_low = y_beam-0.3;
      y_up = y_beam+0.7;
      pT_nbins = 3;
      pT_low = 0.2;
      pT_up = 1.4;
    }
  } else if(std::fabs(y_beam - 0.985344)<0.01) {
    C_binranges = {0, 10, 20, 30, 40, 70};
    if(pdg == 3122) {
      mass_nbins = 1200;
      mass_low = 1.045;
      mass_up = 1.185;
      y_nbins = 7;
      y_low = y_beam-0.3;
      y_up = y_beam+1.1;
      pT_nbins = 3;
      pT_low = 0;
      pT_up = 1.2;
    } else if(pdg == 310) {
      mass_nbins = 1200;
      mass_low = 0.277;
      mass_up = 0.717;
      y_nbins = 7;
      y_low = y_beam-0.1;
      y_up = y_beam+1.3;
      pT_nbins = 3;
      pT_low = 0;
      pT_up = 0.9;
    } else if(pdg == 3312) {
      mass_nbins = 1200;
      mass_low = 1.251;
      mass_up = 1.391;
      y_nbins = 5;
      y_low = y_beam-0.1;
      y_up = y_beam+0.9;
      pT_nbins = 3;
      pT_low = 0;
      pT_up = 0.9;
    }
  }

  pT_nbins *= 2;

  const int C_nbins = C_binranges.size()-1;
  
  Qn::Axis<double> C_axis("centrality", C_binranges);
  Qn::Axis<double> pT_axis("pT", pT_nbins, pT_low, pT_up);
//   Qn::Axis<double> pT_axis("pT", pT_binranges);
  Qn::Axis<double> y_axis("y", y_nbins, y_low, y_up);
  
  Qn::DataContainer<TH1F, Qn::Axis<double>> dcmass_all, dcmass_sgnl, dcmass_bckgr;
  dcmass_all.AddAxes({C_axis, pT_axis, y_axis});
  dcmass_sgnl.AddAxes({C_axis, pT_axis, y_axis});
  dcmass_bckgr.AddAxes({C_axis, pT_axis, y_axis});

  for(int i=0; i<dcmass_all.size(); i++) {
    dcmass_all[i].SetName("hmass");
    dcmass_all[i].SetTitle("hMass");
    dcmass_all[i].SetBins(mass_nbins, mass_low, mass_up);
    dcmass_all[i].GetXaxis()->SetTitle("m_{inv}, GeV");
    dcmass_sgnl[i].SetName("hmass");
    dcmass_sgnl[i].SetTitle("hMass");
    dcmass_sgnl[i].SetBins(mass_nbins, mass_low, mass_up);
    dcmass_sgnl[i].GetXaxis()->SetTitle("m_{inv}, GeV");
    dcmass_bckgr[i].SetName("hmass");
    dcmass_bckgr[i].SetTitle("hMass");
    dcmass_bckgr[i].SetBins(mass_nbins, mass_low, mass_up);
    dcmass_bckgr[i].GetXaxis()->SetTitle("m_{inv}, GeV");
  }
  
  TFile* fileOut = TFile::Open("massDC.root", "recreate");
  
  const int n_entries = treeIn->GetEntries();
  for(int i_event=0; i_event<n_entries; ++i_event)
  {
    if(i_event%100 == 0)
      std::cout << i_event << "\n";
    treeIn->GetEntry(i_event);
    const float centrality = event_header->GetField<float>(centrality_id_);
    for(const auto& reco_track : *(reco_tracks->GetChannels()) )
    {     
      const int pid = reco_track.GetPid();
      if(pid != pdg) continue;
      
      const float pT_reco = reco_track.GetPt();
      const float y_reco = reco_track.GetRapidity();

      const int generation = reco_track.GetField<int>(is_signal_id_);
      
      const int linearindex = dcmass_all.FindBin(std::vector<double>({centrality, pT_reco, y_reco}));
      if(linearindex<0) continue;
      dcmass_all[linearindex].Fill(reco_track.GetMass());
      if(generation == 0) dcmass_bckgr[linearindex].Fill(reco_track.GetMass());
      else                dcmass_sgnl[linearindex].Fill(reco_track.GetMass());
    }    
  }  
      
  dcmass_all.Write("dcmass_all");
  dcmass_bckgr.Write("dcmass_bckgr");
  dcmass_sgnl.Write("dcmass_sgnl");

  fileOut -> Close();
}
