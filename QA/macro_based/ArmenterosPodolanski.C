// v0_ap_sim.C
// A ROOT macro to simulate K_S^0 and Lambda decays and plot Armenteros-Podolanski (AP) variables


void ArmenterosPodolanski(int nEvents = 10000000) {
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
    TH2D* hAP = new TH2D("hAP", "Armenteros-Podolanski; #alpha; p_{T} (GeV/c)",
                        200, -1, 1, 200, 0, 0.3);

    auto doDecay = [&](double M, double m1, double m2, TH2D* h) {
        for (int i = 0; i < nEvents; ++i) {
            // Parent at high momentum along z
            double pParent = rng.Uniform(10, 10); // GeV/c
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

            double alpha = (pL1 - pL2) / (pL1 + pL2);
            double pT = pT1;

            h->Fill(alpha, pT);
        }
    };

//     doDecay(mK0, mpi, mpi, hAP);
//     doDecay(mL, mp, mpi, hAP);
//     doDecay(mL, mpi, mp, hAP);
    doDecay(mOmega, mL, mKpm, hAP);
    doDecay(mOmega, mKpm, mL, hAP);
    doDecay(mXi, mL, mpi, hAP);
    doDecay(mXi, mpi, mL, hAP);

    TCanvas* c = new TCanvas("c", "Armenteros-Podolanski Simulation", 800, 600);
    hAP->Draw("COLZ");
    c->SetLogz();

    c->Print("AP.pdf", "pdf");
}
