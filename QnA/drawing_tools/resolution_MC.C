void resolution_MC() {
  gROOT->Macro( "/home/oleksii/cbmdir/flow_drawing_tools/example/style.cc" );

//   std::string evegen = "dcmqgsm"; std::string pbeam = "12"; std::string cuts = "lc1";
  std::string evegen = "dcmqgsm"; std::string pbeam = "3.3"; std::string cuts = "oc1";
//   std::string evegen = "urqmd";  std::string pbeam = "12"; std::string cuts = "lc1";

//   bool is_write_rootfile = false;
  bool is_write_rootfile = true;
    
//   std::string fileName = "/home/oleksii/cbmdir/working/qna/simtracksflow/" + evegen + "/" + pbeam + "agev/v1andR1.stf." + evegen + "." + pbeam + "agev.root";
  std::string fileName = "/home/oleksii/cbmdir/working/qna/aXmass/cl." + evegen + "." + pbeam + "agev." + cuts + ".3122.root";

  std::vector<std::string> correls{"psd1", "psd2", "psd3"};
//   std::vector<std::string> correls{"etacut_1_charged", "etacut_2_charged", "etacut_3_charged"};
//   std::vector<std::string> correls{"etacut_1_all", "etacut_2_all", "etacut_3_all"};

  std::string step;
  bool average_comp;

//   std::vector<std::string> components{"x1x1", "y1y1"}; std::string L_or_P = "L"; std::string same_or_cross = "res.mc";
  std::vector<std::string> components{"cos1x1", "sin1y1"}; std::string L_or_P = "L"; std::string same_or_cross = "res.mc";
//   std::vector<std::string> components{"x1y1", "y1x1"}; std::string L_or_P = "P"; std::string same_or_cross = "res_cross";

  std::string fileOutName;
  if(correls.at(0)[0] == 'p' && (components.at(0) == "x1x1" || components.at(0) == "cos1x1")) { step = "_RECENTERED"; fileOutName = "res.psd"; average_comp = false; }
  if(correls.at(0)[0] == 'p' && components.at(0) == "x1y1") { step = "_RECENTERED"; fileOutName = "res_cross.psd"; average_comp = false; }
  if(correls.at(0)[0] == 'e' && components.at(0) == "x1x1") { step = "_PLAIN"; fileOutName = "res.etacut"; average_comp = true; }
  if(correls.at(0)[0] == 'e' && components.at(0) == "x1y1") { step = "_PLAIN"; fileOutName = "res_cross.etacut"; average_comp = false; }

  MultiCorrelation multicor_mc;
  multicor_mc.SetIsFillSysErrors(false);
  if(!average_comp) {
    multicor_mc.SetPalette( {kBlue, kBlue, kRed, kRed, kGreen+2, kGreen+2} );
    if(components.at(0) == "x1x1" || components.at(0) == "cos1x1") multicor_mc.SetMarkers( {-1, -2, -1, -2, -1, -2} );
    if(components.at(0) == "x1y1") multicor_mc.SetMarkers( {kFullSquare, kOpenSquare, kFullSquare, kOpenSquare, kFullSquare, kOpenSquare} );
  } else {
    multicor_mc.SetPalette( {kBlue, kRed, kGreen+2} );
    if(components.at(0) == "x1x1" || components.at(0) == "cos1x1") multicor_mc.SetMarkers( {-1, -1, -1} );
  }

  for(auto& corr : correls) {
    if(!average_comp) {
      for(auto& comp : components) {
        multicor_mc.AddCorrelation(fileName, {"QPsi/" + corr + step + ".Q_psi_PLAIN." + comp}, "mc_" + corr + step);
      }
    } else {
      multicor_mc.AddCorrelation(fileName, {"QPsi/" + corr + step + ".Q_psi_PLAIN." + components.at(0),
                                            "QPsi/" + corr + step + ".Q_psi_PLAIN." + components.at(1)}, "mc_" + corr + step);
    }
  }

  multicor_mc.Scale(2);
  multicor_mc.SlightShiftXAxis(0.);
  if(components.at(0) == "x1y1") multicor_mc.SlightShiftXAxis(1);
            
  HeapPicture pic("picture", {1000, 1000});

  const float text_size = 20;

  if(evegen == "dcmqgsm") {
    if(pbeam == "12") pic.AddText("5M Au+Au", text_size, {0.04, 0.96});
    else              pic.AddText("5.2M Au+Au", text_size, {0.04, 0.96});
    pic.AddText("DCM-QGSM-SMM", text_size, {0.04, 0.92});
  }
  if(evegen == "urqmd") {
    pic.AddText("2M Au+Au", text_size, {0.04, 0.96});
    pic.AddText("UrQMD", text_size, {0.04, 0.92});
  }
  pic.AddText(pbeam + "A GeV/c", text_size, {0.04, 0.88});
  pic.AddText("MC: R^{A}_{x} = 2#LTQ^{A}_{x}Q^{#Psi}_{x}#GT",text_size, {0.04, 0.84});

  for(auto& tx : pic.GetTexts()) {
    tx->SetTextFont(63);
  }

  auto leg1 = new TLegend();
  leg1->SetBorderSize(1);

  for(auto& obj : multicor_mc.GetCorrelations()) {
    obj->SetCalculateSystematicsFromVariation(false);
    pic.AddDrawable(obj);
  }

  TGraph* grx = new TGraph();
  grx->SetMarkerStyle(kFullSquare);
  grx->SetLineStyle(1);
  grx->SetMarkerColor(kBlack);
  grx->SetLineColor(kBlack);
  TGraph* gry = new TGraph();
  gry->SetMarkerStyle(kOpenSquare);
  gry->SetLineStyle(2);
  gry->SetMarkerColor(kBlack);
  gry->SetLineColor(kBlack);
  
  if(!average_comp) {
    leg1->AddEntry(grx, components.at(0).c_str(), L_or_P.c_str());
    leg1->AddEntry(gry, components.at(1).c_str(), L_or_P.c_str());
    leg1->AddEntry(multicor_mc.GetCorrelations().at(0)->GetPoints(), (correls.at(0) + step).c_str(), L_or_P.c_str());
    leg1->AddEntry(multicor_mc.GetCorrelations().at(2)->GetPoints(), (correls.at(1) + step).c_str(), L_or_P.c_str());
    leg1->AddEntry(multicor_mc.GetCorrelations().at(4)->GetPoints(), (correls.at(2) + step).c_str(), L_or_P.c_str());
  } else {
    leg1->AddEntry(grx, "x&y ave", L_or_P.c_str());
    leg1->AddEntry(multicor_mc.GetCorrelations().at(0)->GetPoints(), (correls.at(0) + step).c_str(), L_or_P.c_str());
    leg1->AddEntry(multicor_mc.GetCorrelations().at(1)->GetPoints(), (correls.at(1) + step).c_str(), L_or_P.c_str());
    leg1->AddEntry(multicor_mc.GetCorrelations().at(2)->GetPoints(), (correls.at(2) + step).c_str(), L_or_P.c_str());
  }
  
  pic.SetAxisTitles( {"Centrality, %", "R_{1}"} );
  pic.CustomizeXRange();
  pic.CustomizeYRange();  
  pic.AddLegend(leg1);
//   pic.CustomizeLegend(leg1);
//   pic.SetGridX();
//   pic.SetGridY();
  pic.Draw();
//   pic.Save("fileOut", "png");
  
  if(is_write_rootfile) {
    TFile* fileOut = TFile::Open((fileOutName + ".root").c_str(), "recreate");
    fileOut->cd();
    pic.GetCanvas()->Write();
//     pic.Dump();
    pic.Write("heap_picture");
    fileOut->Close();
  }
    
  pic.GetCanvas()->Print((fileOutName + ".pdf").c_str(), "pdf");
}
