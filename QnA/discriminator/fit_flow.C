#include "Helper.hpp"

void fit_flow() {
  gErrorIgnoreLevel = 6001;

  gROOT->Macro("/home/oleksii/cbmdir/flow_drawing_tools/example/style_1.cc");

//   std::string evegen = "dcmqgsm"; std::string pbeam = "12";
  std::string evegen = "urqmd"; std::string pbeam = "12";

  const float mu = 1.115683; std::string particle = "lambda"; std::string pdg = "3122";

//   const float mu = 0.497611; std::string particle = "kshort"; std::string pdg = "310";

//   std::string is_fine_pt = "";
  std::string is_fine_pt = "_finept";

  std::string v1filename = "/home/oleksii/cbmdir/working/qna/aXmass/vR." + evegen + "." + pbeam + "agev.lc1." + pdg + is_fine_pt + ".root";
  std::string shapefilename = "/home/oleksii/cbmdir/working/qna/aXmass/shapes/shape." + evegen + "." + pbeam + "agev.lc1." + pdg + is_fine_pt + ".root";

  TFile* shapefile = TFile::Open(shapefilename.c_str(), "read");
  Qn::DataContainer<Qn::ShapeContainer, Qn::Axis<double>>* shcntr = (Qn::DataContainer<Qn::ShapeContainer, Qn::Axis<double>>*) shapefile->Get(/*"ReFit"*/"PreFit" );

  TFile* v1file = TFile::Open(v1filename.c_str(), "read");

  std::string invmassaxis = "ReconstructedParticles_mass";

  std::vector<std::string> components = {"x1x1", "y1y1"};

  TFile* fileOut = TFile::Open("out.fitter.root", "recreate");

  std::vector<std::pair<std::string, std::string>> inputs {
    {"uPsi", ""},
    {"uQ_R1_MC", "_res_MC"},
    {"uQ_R1_sub3", "_res_sub3"},
    {"uQ_R1_sub4", "_res_sub4_sts_pipos"}
  };

  for(auto& ip : inputs) {
    std::vector<std::string> subevents;
    if(ip.second == "") subevents = {"Q_psi"};
    else                subevents = {"psd1", "psd2", "psd3"};

    for(auto& se : subevents) {
      Qn::DataContainerStatDiscriminator dc_entries_sgnl;
      dc_entries_sgnl.AddAxes(shcntr->GetAxes());

      Qn::DataContainerStatDiscriminator dc_entries_bckgr;
      dc_entries_bckgr.AddAxes(shcntr->GetAxes());

      bool common_dc_for_all_components{false};

      for (auto& co : components) {

        std::string corr_name = "v1/" + ip.first + "/v1.u_rec." + se + ip.second + "." + co;
        std::cout << corr_name << "\n";
        Qn::DataContainer<Qn::StatCalculate, Qn::Axis<double>> lambda_psi_pre =
        Qn::DataContainerStatCalculate(*(Qn::DataContainer<Qn::StatCalculate, Qn::Axis<double>>*) v1file->Get(corr_name.c_str()));

    //     auto lambda_psi = lambda_psi_pre.Rebin({"AnaEventHeader_centrality_tracks", {0, 10, 20, 40, 70}});
        auto lambda_psi = lambda_psi_pre;

        const double invmass_lo = lambda_psi.GetAxis(invmassaxis.c_str()).GetFirstBinEdge();
        const double invmass_hi = lambda_psi.GetAxis(invmassaxis.c_str()).GetLastBinEdge();

        Qn::DataContainer<Qn::StatCalculate, Qn::Axis<double>> lambda_psi_rebinned = lambda_psi.Rebin({invmassaxis, 1, invmass_lo, invmass_hi});
        Qn::DataContainer<Qn::StatCalculate, Qn::Axis<double>> lambda_psi_reduced = lambda_psi_rebinned.Select({invmassaxis, 1, invmass_lo, invmass_hi});

        Qn::DataContainerStatDiscriminator dc_signal;
        dc_signal.AddAxes(shcntr->GetAxes());

        Qn::DataContainerStatDiscriminator dc_bckgr_0;
        dc_bckgr_0.AddAxes(shcntr->GetAxes());

        Qn::DataContainerStatDiscriminator dc_bckgr_1;
        dc_bckgr_1.AddAxes(shcntr->GetAxes());

        Qn::DataContainerStatDiscriminator dc_fit_chi2ndf;
        dc_fit_chi2ndf.AddAxes(shcntr->GetAxes());

        GraphExtractor gex;
        gex.SetDataContainer(&lambda_psi);
        gex.SetSelectAxis(invmassaxis.c_str());

        CD(fileOut, ("Fits/" + ip.first + "/" + se + ip.second).c_str());

        const int Nsamples = lambda_psi.At(0).GetSampleMeans().size();

        for (int i = 0; i < shcntr->size(); i++) {
          gex.ReduceDataContainerToBin(shcntr->GetIndex(i));
          TGraphErrors* gr = gex.GetGraph();
          std::vector<double> samples_weights = gex.GetSamplesWeights();
          std::vector<TGraphErrors*> grs = gex.GetSamplesGraphs();
          std::vector indices = shcntr->GetIndex(i);
          std::string binname = "C" + StringBinNumber(indices.at(0) + 1) + "_pT" + StringBinNumber(indices.at(1) + 1) + "_y" + StringBinNumber(indices.at(2) + 1) + "." + co;
          const float C_lo = shcntr->GetAxis("centrality").GetLowerBinEdge(indices.at(0));
          const float C_hi = shcntr->GetAxis("centrality").GetUpperBinEdge(indices.at(0));
          const float pT_lo = shcntr->GetAxis("pT").GetLowerBinEdge(indices.at(1));
          const float pT_hi = shcntr->GetAxis("pT").GetUpperBinEdge(indices.at(1));
          const float y_lo = shcntr->GetAxis("y").GetLowerBinEdge(indices.at(2));
          const float y_hi = shcntr->GetAxis("y").GetUpperBinEdge(indices.at(2));

          std::cout << "DataContainer cell # " << i << "\n";
          std::cout << binname << "\n";
          std::cout << ("C: " + to_string_with_precision(C_lo, 2) + " - " + to_string_with_precision(C_hi, 2) + " %").c_str() << "\n";
          std::cout << ("p_{T}: " + to_string_with_precision(pT_lo, 2) + " - " + to_string_with_precision(pT_hi, 2) + " GeV/c").c_str() << "\n";
          std::cout << ("y_{LAB}: " + to_string_with_precision(y_lo, 2) + " - " + to_string_with_precision(y_hi, 2)).c_str() << "\n";

          const float sumweights = lambda_psi_reduced[i].SumWeights();
          const float sgnl_integral = shcntr->At(i).GetSignalIntegral(invmass_lo, invmass_hi);
          const float bckgr_integral = shcntr->At(i).GetBackgroundIntegral(invmass_lo, invmass_hi);

          const float sgnl_weight = sumweights * sgnl_integral / (sgnl_integral + bckgr_integral);
          const float bckgr_weight = sumweights * bckgr_integral / (sgnl_integral + bckgr_integral);

          Fitter fitter;
          fitter.SetMu(mu);
          fitter.SetShape(&shcntr->At(i));
          fitter.SetGraphToFit(gr);
          fitter.SetBsGraphsToFit(grs);
          fitter.Fit();
    //       fitter.Print();

          TCanvas c1("", "", 1500, 900);
          c1.cd();
          gr->SetName(binname.c_str());
          gr->GetXaxis()->SetTitle("m_{inv}, GeV");
          gr->GetYaxis()->SetTitle("v_{1}");
          gr->SetTitle(binname.c_str());
          gr->Draw("AP");
          fitter.GetBckgrGraph()->SetFillStyle(3001);
          fitter.GetBckgrGraph()->SetFillColor(kAzure + 8);
          fitter.GetBckgrGraph()->SetLineColor(kBlue);
          fitter.GetBckgrGraph()->Draw("l e3 same");
          fitter.GetGraphFit()->SetFillStyle(3001);
          fitter.GetGraphFit()->SetFillColor(kRed - 4);
          fitter.GetGraphFit()->SetLineColor(kRed);
          fitter.GetGraphFit()->Draw("l e3 same");

          TPaveText binedges(0.10, 0.78, 0.25, 0.92, "brNDC");
          binedges.AddText(("C: " + to_string_with_precision(C_lo, 2) + " - " + to_string_with_precision(C_hi, 2) + " %").c_str());
          binedges.AddText(("p_{T}: " + to_string_with_precision(pT_lo, 2) + " - " + to_string_with_precision(pT_hi, 2) + " GeV/c").c_str());
          binedges.AddText(("y_{LAB}: " + to_string_with_precision(y_lo, 2) + " - " + to_string_with_precision(y_hi, 2)).c_str());
          binedges.SetFillColor(0);
          binedges.SetTextSize(0.03);
          binedges.SetTextFont(22);
          binedges.Draw("same");

          TLegend legend(0.12, 0.62, 0.27, 0.73);
          legend.SetBorderSize(0);
          legend.AddEntry(fitter.GetGraphFit(), "All Fit", "L");
          legend.AddEntry(fitter.GetBckgrGraph(), "BG Fit", "L");
          legend.SetTextSize(0.03);
          legend.SetTextFont(22);
          legend.Draw("same");

          const float par_sgnl = fitter.GetFitParameters().at(0);
          const float parerr_sgnl = fitter.GetFitErrors().at(0);
          const int npar_bckgr = 2;
          std::vector<float> par_bckgr;
          std::vector<float> parerr_bckgr;
          for (int i = 0; i < npar_bckgr; i++) {
            par_bckgr.push_back(fitter.GetFitParameters().at(i + 1));
            parerr_bckgr.push_back(fitter.GetFitErrors().at(i + 1));
          }
          TPaveText ptpar(0.76, 0.77, 0.91, 0.92, "brNDC");
          ptpar.AddText(("v_{S} = " + to_string_with_precision(par_sgnl, 3) + " #pm " + to_string_with_precision(parerr_sgnl, 3)).c_str());
          ptpar.AddText(("v_{BG} |_{m_{inv}=#mu} = " + to_string_with_precision(par_bckgr.at(0), 3) + " #pm " + to_string_with_precision(parerr_bckgr.at(0), 3)).c_str());
          ptpar.AddText(("dv_{BG}/dm_{inv} |_{m_{inv}=#mu} = " + to_string_with_precision(par_bckgr.at(1), 3) + " #pm " + to_string_with_precision(parerr_bckgr.at(1), 3)).c_str());
          ptpar.SetFillColor(0);
          ptpar.SetTextSize(0.03);
          ptpar.SetTextFont(22);
          ptpar.Draw("same");

          const float chi2_fit = fitter.GetFitChi2Ndf();
          TPaveText ptchi2(0.30, 0.82, 0.45, 0.92, "brNDC");
          ptchi2.AddText(("#chi^{2}_{ndf} Fit = " + to_string_with_precision(chi2_fit, 2)).c_str());
          ptchi2.SetFillColor(0);
          ptchi2.SetTextSize(0.03);
          ptchi2.SetTextFont(22);
          ptchi2.Draw("same");

          TPaveText ptyield(0.12, 0.10, 0.27, 0.20, "brNDC");
          ptyield.AddText(("N_{sgnl} = " + to_string_with_precision(sgnl_weight, 0)).c_str());
          ptyield.AddText(("N_{bckgr} = " + to_string_with_precision(bckgr_weight, 0)).c_str());
          ptyield.SetFillColor(0);
          ptyield.SetTextSize(0.03);
          ptyield.SetTextFont(22);
          ptyield.Draw("same");

          c1.Write(gr->GetName());

          dc_signal[i].SetValue(par_sgnl);
          dc_signal[i].SetError(parerr_sgnl);
          dc_signal[i].SetWeight(sgnl_weight);

          dc_bckgr_0[i].SetValue(par_bckgr.at(0));
          dc_bckgr_0[i].SetError(parerr_bckgr.at(0));
          dc_bckgr_0[i].SetWeight(bckgr_weight);

          dc_bckgr_1[i].SetValue(par_bckgr.at(1));
          dc_bckgr_1[i].SetError(parerr_bckgr.at(1));
          dc_bckgr_1[i].SetWeight(bckgr_weight);

          dc_fit_chi2ndf[i].SetVEW(chi2_fit);

          for(int isample=0; isample<Nsamples; isample++) {
            dc_signal[i].AddSampleMean(fitter.GetBsFitParameter(isample, 0));
            dc_bckgr_0[i].AddSampleMean(fitter.GetBsFitParameter(isample, 1));
            dc_bckgr_1[i].AddSampleMean(fitter.GetBsFitParameter(isample, 2));
            dc_signal[i].AddSampleWeight(samples_weights.at(isample));        //TODO check if it is correct
            dc_bckgr_0[i].AddSampleWeight(samples_weights.at(isample));
            dc_bckgr_1[i].AddSampleWeight(samples_weights.at(isample));
          }

          std::cout << "\n\n";

          if (common_dc_for_all_components) continue;

          dc_entries_sgnl[i].SetValue(sgnl_weight);
          dc_entries_sgnl[i].SetError(0);
          dc_entries_sgnl[i].SetWeight(sgnl_weight);

          dc_entries_bckgr[i].SetValue(bckgr_weight);
          dc_entries_bckgr[i].SetError(0);
          dc_entries_bckgr[i].SetWeight(bckgr_weight);
        }

        CD(fileOut, ("Pars/" + ip.first + "/" + se + ip.second).c_str());
        dc_signal.Write(("signal." + co).c_str());
        dc_bckgr_0.Write(("bckgr_0." + co).c_str());
        dc_bckgr_1.Write(("bckgr_1." + co).c_str());

        CD(fileOut, ("Chi2s/" + ip.first + "/" + se + ip.second).c_str());
        dc_fit_chi2ndf.Write(("fit_chi2ndf." + co).c_str());

        common_dc_for_all_components = true;
      } // components

      CD(fileOut, ("Entries/" + ip.first + "/" + se + ip.second).c_str());
      dc_entries_sgnl.Write("entries_sgnl");
      dc_entries_bckgr.Write("entries_bckgr");


    } // subevents
  } // inputs

  fileOut->Close();
}
