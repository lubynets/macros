// v0_ap_sim.C
// A ROOT macro to simulate K_S^0 and Lambda decays and plot Armenteros-Podolanski (AP) variables

double evalAlpha0(double M, double m1, double m2);
double evalQStar(double M, double m1, double m2);
double evalRalpha(double M, double m1, double m2, double P);
bool checkCascade(double aplpha, double qt);

void ArmenterosPodolanskiToyMC(int nEvents = 10000000) {
    TRandom3 rng(0);

    // Masses in GeV/c^2
    const double mK0 = 0.497611;      // K_S^0
    const double mpi = 0.139570;      // charged pion
    const double mL = 1.115683;       // Lambda
    const double mp = 0.938272;       // proton
    const double mOmega = 1.67245;    // omega
    const double mKpm = 0.493677;     // K^+-
    const double mXi = 1.32171;       // Xi^-

    // Histograms: qT vs alpha
    TH2D* hAP = new TH2D("hAP", "Armenteros-Podolanski; #alpha; p_{T} (GeV/c)", 200, -1, 1, 200, 0, 0.3);

    TH2D* hPBachAlpha = new TH2D("hPBachAlpha", "", 200, 0.1, 10, 200, -1, 1);
    TH2D* hPBachQt = new TH2D("hPBachQt", "", 200, 0.1, 10, 200, 0, 0.3);
    hPBachAlpha->GetXaxis()->SetTitle("#it{p}_{bach} (GeV/#it{c})");
    hPBachAlpha->GetYaxis()->SetTitle("#alpha");
    hPBachQt->GetXaxis()->SetTitle("#it{p}_{bach} (GeV/#it{c})");
    hPBachQt->GetYaxis()->SetTitle("#it{q}_{T}");


    auto doDecay = [&](double M, double m1, double m2, TH2D* hAP, TH2D* hPBachAlpha=nullptr, TH2D* hPBachQt=nullptr) {
      std::cout << "alpha0 = " << evalAlpha0(M, m1, m2) << "\n";
      std::cout << "QStar = " << evalQStar(M, m1, m2) << "\n";
      std::cout << "evalRalpha(P=1) = " << evalRalpha(M, m1, m2, 1) << "\n";
      std::cout << "evalRalpha(P=5) = " << evalRalpha(M, m1, m2, 5) << "\n";
      std::cout << "evalRalpha(P=10) = " << evalRalpha(M, m1, m2, 10) << "\n";
      for (int i = 0; i < nEvents; ++i) {
        // Parent at high momentum along z
        double pParent = rng.Uniform(1, 10); // GeV/c
        TLorentzVector P(0, 0, pParent, sqrt(pParent*pParent + M*M));

        // Daughter momenta in parent rest frame (two-body decay)
        double pStar = TMath::Sqrt((M*M - (m1 + m2)*(m1 + m2)) * (M*M - (m1 - m2)*(m1 - m2))) / (2*M);
        double costh = rng.Uniform(-1, 1);
        double phi = rng.Uniform(0, 2*TMath::Pi());
        TLorentzVector p1_rf(
            pStar * TMath::Sin(acos(costh)) * TMath::Cos(phi),
            pStar * TMath::Sin(acos(costh)) * TMath::Sin(phi),
            pStar * costh,
            sqrt(pStar*pStar + m1*m1)
        );
        TLorentzVector p2_rf(-p1_rf.Vect(), sqrt(pStar*pStar + m2*m2));

        // Boost daughters to lab frame
        TVector3 beta = P.BoostVector();
        p1_rf.Boost(beta);
        p2_rf.Boost(beta);

        double pL1 = p1_rf.Vect().Dot(P.Vect().Unit());
        double pL2 = p2_rf.Vect().Dot(P.Vect().Unit());
        double pT1 = p1_rf.Vect().Perp(P.Vect().Unit());
        // pT2 is the same by momentum conservation
        const double pBach = p2_rf.Vect().Mag();

        double alpha = (pL1 - pL2) / (pL1 + pL2);
//         if(alpha < 0.05 || alpha > 0.5) continue;
        double pT = pT1;

        hAP->Fill(alpha, pT);
        if(hPBachAlpha != nullptr && hPBachQt != nullptr) {
          hPBachAlpha->Fill(pBach, alpha);
          hPBachQt->Fill(pBach, pT);
        }
      }
    };

//     doDecay(mK0, mpi, mpi, hAP);
//     doDecay(mL, mp, mpi, hAP);
//     doDecay(mL, mpi, mp, hAP);
    doDecay(mOmega, mL, mKpm, hAP, hPBachAlpha, hPBachQt);
//     doDecay(mOmega, mKpm, mL, hAP);
//     doDecay(mXi, mL, mpi, hAP);
//     doDecay(mXi, mpi, mL, hAP);

    TFile* fileOut = TFile::Open("fileOut.root", "recreate");
    hAP->Write();
    hPBachAlpha->Write();
    hPBachQt->Write();
    fileOut->Close();

    TCanvas* c = new TCanvas("c", "Armenteros-Podolanski Simulation", 800, 600);
    hAP->Draw("COLZ");
    c->SetLogz();

    c->Print("AP.pdf", "pdf");
}

double evalAlpha0(double M, double m1, double m2) {
  return std::fabs(m1*m1 - m2*m2) / M / M;
}

double evalQStar(double M, double m1, double m2) {
  return std::sqrt(M*M*M*M + m1*m1*m1*m1 + m2*m2*m2*m2 + - 2.*M*M*(m1*m1 + m2*m2) - 2.*m1*m1*m2*m2) / 2. / M;
}

double evalRalpha(double M, double m1, double m2, double P) {
  const double b = P / std::sqrt(P*P + M*M);
  return 2*evalQStar(M, m1, m2) / b / M;
}
