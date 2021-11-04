#include "Particle.h"

#include "Parameters.h"
#include "ParticleType.h"
#include "ResonanceType.h"

#include <cmath>
#include <iostream>
#include <string>
#include <TH1C.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TLatex.h>

void Main() {
    // Initialization of particle types
    Particle::AddParticleType("π+", 0.13957, +1);
    Particle::AddParticleType("π-", 0.13957, -1);
    Particle::AddParticleType("K+", 0.49367, +1);
    Particle::AddParticleType("K-", 0.49367, -1);
    Particle::AddParticleType("p+", 0.93827, +1);
    Particle::AddParticleType("p-", 0.93827, -1);
    Particle::AddParticleType("K*", 0.89166, 0, 0.050);

    // Histograms definitions
    TList *histograms = new TList();

    TH1F *particleTypesH = new TH1F("particleTypesH", "Particle Types", N_PARTICLE_TYPES, 0, N_PARTICLE_TYPES);
    particleTypesH->SetFillColor(kBlue);
    particleTypesH->GetXaxis()->SetBinLabel(PION_PLUS_BIN, "#pi+");
    particleTypesH->GetXaxis()->SetBinLabel(PION_MINUS_BIN, "#pi-");
    particleTypesH->GetXaxis()->SetBinLabel(KAON_PLUS_BIN, "K+");
    particleTypesH->GetXaxis()->SetBinLabel(KAON_MINUS_BIN, "K-");
    particleTypesH->GetXaxis()->SetBinLabel(PROTON_PLUS_BIN, "p+");
    particleTypesH->GetXaxis()->SetBinLabel(PROTON_MINUS_BIN, "p+");
    particleTypesH->GetXaxis()->SetBinLabel(KAON_STAR_BIN, "K*");
    particleTypesH->SetStats(false);
    histograms->Add(particleTypesH);

    TH1D *azimutAngleH = new TH1D("azimutAngleH", "Azimut Angle Distribution", N_BINS, 0, 2 * M_PI);
    histograms->Add(azimutAngleH);

    TH1D *polarAngleH = new TH1D("polarAngleH", "Polar Angle Distribution", N_BINS, 0, M_PI);
    histograms->Add(polarAngleH);

    TH1D *momentumH = new TH1D("momentumH", "Momentum Distribution", N_BINS, 0, MAX_MOMENTUM);
    histograms->Add(momentumH);

    TH1D *transverseMomentumH = new TH1D("transverseMomentumH", "Transverse Momentum Distribution", N_BINS, 0, MAX_MOMENTUM);
    histograms->Add(transverseMomentumH);

    TH1D *particleEnergyH = new TH1D("particleEnergyH", "Particle Energy Distribution", N_BINS, 0, MAX_ENERGY);
    histograms->Add(particleEnergyH);

    TH1D *invMassH = new TH1D("invMassH", "Invariant Mass Distribution", N_BINS, 0, MAX_INVARIANT_MASS);
    invMassH->Sumw2();
    histograms->Add(invMassH);

    TH1D *discordantInvMassH = new TH1D("discordantInvMassH", "Discordant Invariant Mass Distribution", N_BINS, 0, MAX_INVARIANT_MASS);
    discordantInvMassH->Sumw2();
    histograms->Add(discordantInvMassH);

    TH1D *concordantInvMassH = new TH1D("concordantInvMassH", "Concordant Invariant Mass Distribution", N_BINS, 0, MAX_INVARIANT_MASS);
    concordantInvMassH->Sumw2();
    histograms->Add(concordantInvMassH);

    TH1D *discordantPionKaonInvMassH = new TH1D("discordantPionKaonInvMassH", "Discordant Pion/Kaon Invariant Mass Distribution", N_BINS, 0, MAX_INVARIANT_MASS);
    discordantPionKaonInvMassH->Sumw2();
    histograms->Add(discordantPionKaonInvMassH);

    TH1D *concordantPionKaonInvMassH = new TH1D("concordantPionKaonInvMassH", "Concordant Pion/Kaon Invariant Mass Distribution", N_BINS, 0, MAX_INVARIANT_MASS);
    concordantPionKaonInvMassH->Sumw2();
    histograms->Add(concordantPionKaonInvMassH);

    TH1D *daughtersInvMassH = new TH1D("daughtersInvMassH", "K* Daughters Invariant Mass Distribution", N_BINS, 0, MAX_INVARIANT_MASS);
    daughtersInvMassH->Sumw2();
    histograms->Add(daughtersInvMassH);

    // Variables definitions
    Particle particles[N_PARTICLE_TYPES + MAX_PRODUCTS];
    double phi, theta, P, rndm;
    double Px, Py, Pz;
    double invMass;

    int nDecayedParticles;  // Counter of decayed particles

    for (int i = 0; i < N_ITERATIONS; i++) {
        nDecayedParticles = 0;
        for (int j = 0; j < N_PARTICLES; j++) {
            // Random generation of momentum
            phi = gRandom->Uniform(0, 2*M_PI);
            theta = gRandom->Uniform(0, M_PI);
            P = gRandom->Exp(AVG_P);

            Px = P * sin(theta) * cos(phi);
            Py = P * sin(theta) * sin(phi);
            Pz = P * cos(theta);
            particles[j].SetP(Px, Py, Pz);

            // Random generation of particle type
            rndm = gRandom->Rndm();
            if (rndm < PION_PLUS_CUMULATIVE) {
                particles[j].SetIndex("π+");
                particleTypesH->Fill(particleTypesH->GetBinCenter(PION_PLUS_BIN));
            } else if (rndm < PION_MINUS_CUMULATIVE) {
                particles[j].SetIndex("π-");
                particleTypesH->Fill(particleTypesH->GetBinCenter(PION_MINUS_BIN));
            } else if (rndm < KAON_PLUS_CUMULATIVE) {
                particles[j].SetIndex("K+");
                particleTypesH->Fill(particleTypesH->GetBinCenter(KAON_PLUS_BIN));
            } else if (rndm < KAON_MINUS_CUMULATIVE) {
                particles[j].SetIndex("K-");
                particleTypesH->Fill(particleTypesH->GetBinCenter(KAON_MINUS_BIN));
            } else if (rndm < PROTON_PLUS_CUMULATIVE) {
                particles[j].SetIndex("p+");
                particleTypesH->Fill(particleTypesH->GetBinCenter(PROTON_PLUS_BIN));
            } else if (rndm < PROTON_MINUS_CUMULATIVE) {
                particles[j].SetIndex("p-");
                particleTypesH->Fill(particleTypesH->GetBinCenter(PROTON_MINUS_BIN));
            } else {
                particles[j].SetIndex("K*");
                particleTypesH->Fill(particleTypesH->GetBinCenter(KAON_STAR_BIN));

                // Decayment of K*
                rndm = gRandom->Rndm();
                if (rndm < 0.5) {
                    particles[N_PARTICLES + 2 * nDecayedParticles].SetIndex("π+");
                    particles[N_PARTICLES + 2 * nDecayedParticles + 1].SetIndex("K-");
                } else {
                    particles[N_PARTICLES + 2 * nDecayedParticles].SetIndex("π-");
                    particles[N_PARTICLES + 2 * nDecayedParticles + 1].SetIndex("K+");
                }
                particles[j].Decay2Body(particles[N_PARTICLES + 2 * nDecayedParticles], particles[N_PARTICLES + 2 * nDecayedParticles + 1]);
                nDecayedParticles++;
            }

            azimutAngleH->Fill(phi);
            polarAngleH->Fill(theta);
            momentumH->Fill(P);
            transverseMomentumH->Fill(P * sin(theta));
            particleEnergyH->Fill(particles[j].TotEnergy());
        }

        // Compute invariant masses and fill histograms
        for (int j = 0; j < N_PARTICLES + 2 * nDecayedParticles; j++) {
            const ParticleType *p1 = particles[j].GetParticleType();

            // Ignore K* particles for invariant mass histograms (charge == 0)
            if (p1->GetCharge() != 0) {
                for (int k = 0; k < j; k++) {
                    const ParticleType *p2 = particles[k].GetParticleType();
                    if (p2->GetCharge() != 0) {
                        invMass = particles[j].InvMass(&particles[k]);
                        invMassH->Fill(invMass);

                        if (p1->GetCharge() != p2->GetCharge())
                            discordantInvMassH->Fill(invMass);
                        else
                            concordantInvMassH->Fill(invMass);

                        if ((p1->GetName() == "π+" && p2->GetName() == "K-") ||
                            (p1->GetName() == "π-" && p2->GetName() == "K+"))
                            discordantPionKaonInvMassH->Fill(invMass);
                        else if ((p1->GetName() == "π+" && p2->GetName() == "K+") ||
                                 (p1->GetName() == "π-" && p2->GetName() == "K-"))
                            concordantPionKaonInvMassH->Fill(invMass);
                    }
                }
            }

            // Fill histogram with invariant masses of K* daughters
            if (j >= N_PARTICLES && j % 2 == 0)
                daughtersInvMassH->Fill(particles[j].InvMass(&particles[j + 1]));
        }
    }

    // Check histograms and their properties
    if (particleTypesH->GetEntries() != N_PARTICLES * N_ITERATIONS)
        std::cout << "Number of entries of Particle Types Histogram is incorrect" << std::endl;
    if (azimutAngleH->GetEntries() != N_PARTICLES * N_ITERATIONS)
        std::cout << "Number of entries of Azimut Angle Histogram is incorrect" << std::endl;
    if (polarAngleH->GetEntries() != N_PARTICLES * N_ITERATIONS)
        std::cout << "Number of entries of Polar Angle Histogram is incorrect" << std::endl;
    if (momentumH->GetEntries() != N_PARTICLES * N_ITERATIONS)
        std::cout << "Number of entries of Momentum Histogram is incorrect" << std::endl;
    if (transverseMomentumH->GetEntries() != N_PARTICLES * N_ITERATIONS)
        std::cout << "Number of entries of Transverse Momentum Histogram is incorrect" << std::endl;
    if (particleEnergyH->GetEntries() != N_PARTICLES * N_ITERATIONS)
        std::cout << "Number of entries of Particle Energy Histogram is incorrect" << std::endl;

    const double nFinalParticlesPerIteration = (N_PARTICLES * N_ITERATIONS + particleTypesH->GetBinContent(KAON_STAR_BIN)) / N_ITERATIONS;
    const double nFinalParticlesPerIterationErr = particleTypesH->GetBinError(KAON_STAR_BIN) / N_ITERATIONS;
    const double nInvMassEntries = (nFinalParticlesPerIteration * (nFinalParticlesPerIteration - 1) / 2) * N_ITERATIONS;
    const double nInvMassEntriesErr = 2 * nFinalParticlesPerIterationErr * nInvMassEntries / nFinalParticlesPerIteration;
    if (invMassH->GetEntries() < nInvMassEntries - ERROR_FACTOR * nInvMassEntriesErr || invMassH->GetEntries() > nInvMassEntries + ERROR_FACTOR * nInvMassEntriesErr)
        std::cout << "Number of entries of Invariant Mass Histogram is incorrect" << std::endl;
    /*
    if (discordantInvMassH->GetEntries() != N_PARTICLES * N_ITERATIONS)
        std::cout << "Number of entries of Particle Types Histogram is incorrect" << std::endl;
    if (concordantInvMassH->GetEntries() != N_PARTICLES * N_ITERATIONS)
        std::cout << "Number of entries of Particle Types Histogram is incorrect" << std::endl;
    if (discordantPionKaonInvMassH->GetEntries() != N_PARTICLES * N_ITERATIONS)
        std::cout << "Number of entries of Particle Types Histogram is incorrect" << std::endl;
    if (concordantPionKaonInvMassH->GetEntries() != N_PARTICLES * N_ITERATIONS)
        std::cout << "Number of entries of Particle Types Histogram is incorrect" << std::endl;
    if (daughtersInvMassH->GetEntries() != N_PARTICLES * N_ITERATIONS)
        std::cout << "Number of entries of Particle Types Histogram is incorrect" << std::endl;
    */

    // Check proportions of particle types
    for (int i = 1; i <= particleTypesH->GetNbinsX(); i++) {
        const double nProb = N_PARTICLES * N_ITERATIONS * PROBABILITIES[i - 1];
        const double nParticles = particleTypesH->GetBinContent(i);
        const double nParticlesErr = particleTypesH->GetBinError(i);
        if (nProb < nParticles - ERROR_FACTOR * nParticlesErr || nProb > nParticles + ERROR_FACTOR * nParticlesErr)
            std::cout << "Number of " << LABELS[i - 1] << " particles is incorrect" << std::endl;
    }

    // Fit histograms
    TFitResultPtr azimutAngleFit = azimutAngleH->Fit("pol0", "S");
    TFitResultPtr polarAngleFit = polarAngleH->Fit("pol0", "S");
    TFitResultPtr momentumFit = momentumH->Fit("expo", "S");

    // Check consistency of invariant mass histograms
    TH1D *pionKaonDiscordantMinusConcordantH = (TH1D*) discordantPionKaonInvMassH->Clone("pionKaonDiscordantMinusConcordantH");
    pionKaonDiscordantMinusConcordantH->Add(concordantPionKaonInvMassH, -1.0);
    histograms->Add(pionKaonDiscordantMinusConcordantH);
    TFitResultPtr pionKaonDiscordantMinusConcordantFit = pionKaonDiscordantMinusConcordantH->Fit("gaus", "S");

    TH1D *discordantMinusConcordantH = (TH1D*) discordantInvMassH->Clone("discordantMinusConcordantH");
    discordantMinusConcordantH->Add(concordantInvMassH, -1.0);
    histograms->Add(discordantMinusConcordantH);
    TFitResultPtr discordantMinusConcordantFit = discordantMinusConcordantH->Fit("gaus", "S");

    // Save histograms in root files
    TList *canvases = new TList();
    TFile *file = new TFile("data.root", "RECREATE");
    for (TObject *&&histogram : *histograms)
        histogram->Write();
    file->Close();

    // Clean up memory
    delete particleTypesH;
    delete azimutAngleH;
    delete polarAngleH;
    delete momentumH;
    delete transverseMomentumH;
    delete particleEnergyH;
    delete invMassH;
    delete discordantInvMassH;
    delete concordantInvMassH;
    delete discordantPionKaonInvMassH;
    delete concordantPionKaonInvMassH;
    delete daughtersInvMassH;
}