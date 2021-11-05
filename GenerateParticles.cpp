#include "Particle.h"

#include "Parameters.h"
#include "ParticleType.h"
#include "ResonanceType.h"

#include <cmath>
#include <iostream>
#include <TH1D.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TROOT.h>

void GenerateParticles() {
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
    
    // Save histograms in root file
    TFile *file = new TFile("histograms.root", "RECREATE");

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

    TH1F *finalParticleTypesH = new TH1F("finalParticleTypesH", "Final Particle Types", N_PARTICLE_TYPES, 0, N_PARTICLE_TYPES);
    finalParticleTypesH->SetFillColor(kBlue);
    finalParticleTypesH->GetXaxis()->SetBinLabel(PION_PLUS_BIN, "#pi+");
    finalParticleTypesH->GetXaxis()->SetBinLabel(PION_MINUS_BIN, "#pi-");
    finalParticleTypesH->GetXaxis()->SetBinLabel(KAON_PLUS_BIN, "K+");
    finalParticleTypesH->GetXaxis()->SetBinLabel(KAON_MINUS_BIN, "K-");
    finalParticleTypesH->GetXaxis()->SetBinLabel(PROTON_PLUS_BIN, "p+");
    finalParticleTypesH->GetXaxis()->SetBinLabel(PROTON_MINUS_BIN, "p+");
    finalParticleTypesH->GetXaxis()->SetBinLabel(KAON_STAR_BIN, "K*");
    finalParticleTypesH->SetStats(false);
    histograms->Add(finalParticleTypesH);

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

    TH1D *invMassH = new TH1D("invMassH", "Invariant Mass Distribution", N_BINS, MIN_INVARIANT_MASS, MAX_INVARIANT_MASS);
    invMassH->Sumw2();
    histograms->Add(invMassH);

    TH1D *discordantInvMassH = new TH1D("discordantInvMassH", "Discordant Invariant Mass Distribution", N_BINS, MIN_INVARIANT_MASS, MAX_INVARIANT_MASS);
    discordantInvMassH->Sumw2();
    histograms->Add(discordantInvMassH);

    TH1D *concordantInvMassH = new TH1D("concordantInvMassH", "Concordant Invariant Mass Distribution", N_BINS, MIN_INVARIANT_MASS, MAX_INVARIANT_MASS);
    concordantInvMassH->Sumw2();
    histograms->Add(concordantInvMassH);

    TH1D *discordantPionKaonInvMassH = new TH1D("discordantPionKaonInvMassH", "Discordant Pion/Kaon Invariant Mass Distribution", N_BINS, MIN_INVARIANT_MASS, MAX_INVARIANT_MASS);
    discordantPionKaonInvMassH->Sumw2();
    histograms->Add(discordantPionKaonInvMassH);

    TH1D *concordantPionKaonInvMassH = new TH1D("concordantPionKaonInvMassH", "Concordant Pion/Kaon Invariant Mass Distribution", N_BINS, MIN_INVARIANT_MASS, MAX_INVARIANT_MASS);
    concordantPionKaonInvMassH->Sumw2();
    histograms->Add(concordantPionKaonInvMassH);

    TH1D *daughtersInvMassH = new TH1D("daughtersInvMassH", "K* Daughters Invariant Mass Distribution", N_BINS, MIN_INVARIANT_MASS, MAX_INVARIANT_MASS);
    daughtersInvMassH->Sumw2();
    histograms->Add(daughtersInvMassH);

    // Variable definitions
    Particle particles[N_PARTICLE_TYPES + MAX_PRODUCTS];
    double phi, theta, P, rndm;
    double Px, Py, Pz;
    double invMass;

    int nDecayedParticles;  // Counter of decayed particles

    // Iterations
    for (int i = 0; i < N_ITERATIONS; i++) {

        // Reset decayed particles counter from previous iterations
        nDecayedParticles = 0;

        // Fill particles array
        for (int j = 0; j < N_PARTICLES_PER_ITERATION; j++) {

            // Random generation of momentum
            phi = gRandom->Uniform(0, 2*M_PI);
            theta = gRandom->Uniform(0, M_PI);
            P = gRandom->Exp(AVG_P);

            Px = P * sin(theta) * cos(phi);
            Py = P * sin(theta) * sin(phi);
            Pz = P * cos(theta);
            particles[j].SetP(Px, Py, Pz);

            // Random generate particle type and fill correspondent histogram
            rndm = gRandom->Rndm();
            if (rndm < PION_PLUS_CUMULATIVE) {
                particles[j].SetIndex("π+");
                particleTypesH->Fill(particleTypesH->GetBinCenter(PION_PLUS_BIN));
                finalParticleTypesH->Fill(finalParticleTypesH->GetBinCenter(PION_PLUS_BIN));
            } else if (rndm < PION_MINUS_CUMULATIVE) {
                particles[j].SetIndex("π-");
                particleTypesH->Fill(particleTypesH->GetBinCenter(PION_MINUS_BIN));
                finalParticleTypesH->Fill(finalParticleTypesH->GetBinCenter(PION_MINUS_BIN));
            } else if (rndm < KAON_PLUS_CUMULATIVE) {
                particles[j].SetIndex("K+");
                particleTypesH->Fill(particleTypesH->GetBinCenter(KAON_PLUS_BIN));
                finalParticleTypesH->Fill(finalParticleTypesH->GetBinCenter(KAON_PLUS_BIN));
            } else if (rndm < KAON_MINUS_CUMULATIVE) {
                particles[j].SetIndex("K-");
                particleTypesH->Fill(particleTypesH->GetBinCenter(KAON_MINUS_BIN));
                finalParticleTypesH->Fill(finalParticleTypesH->GetBinCenter(KAON_MINUS_BIN));
            } else if (rndm < PROTON_PLUS_CUMULATIVE) {
                particles[j].SetIndex("p+");
                particleTypesH->Fill(particleTypesH->GetBinCenter(PROTON_PLUS_BIN));
                finalParticleTypesH->Fill(finalParticleTypesH->GetBinCenter(PROTON_PLUS_BIN));
            } else if (rndm < PROTON_MINUS_CUMULATIVE) {
                particles[j].SetIndex("p-");
                particleTypesH->Fill(particleTypesH->GetBinCenter(PROTON_MINUS_BIN));
                finalParticleTypesH->Fill(finalParticleTypesH->GetBinCenter(PROTON_MINUS_BIN));
            } else {
                particles[j].SetIndex("K*");
                particleTypesH->Fill(particleTypesH->GetBinCenter(KAON_STAR_BIN));
                finalParticleTypesH->Fill(finalParticleTypesH->GetBinCenter(KAON_STAR_BIN));

                // Decayment of K* in random pair (π+, K-) or (π-, K+)
                rndm = gRandom->Rndm();
                if (rndm < 0.5) {
                    particles[N_PARTICLES_PER_ITERATION + 2 * nDecayedParticles].SetIndex("π+");
                    particles[N_PARTICLES_PER_ITERATION + 2 * nDecayedParticles + 1].SetIndex("K-");
                    finalParticleTypesH->Fill(finalParticleTypesH->GetBinCenter(PION_PLUS_BIN));
                    finalParticleTypesH->Fill(finalParticleTypesH->GetBinCenter(KAON_MINUS_BIN));
                } else {
                    particles[N_PARTICLES_PER_ITERATION + 2 * nDecayedParticles].SetIndex("π-");
                    particles[N_PARTICLES_PER_ITERATION + 2 * nDecayedParticles + 1].SetIndex("K+");
                    finalParticleTypesH->Fill(finalParticleTypesH->GetBinCenter(PION_MINUS_BIN));
                    finalParticleTypesH->Fill(finalParticleTypesH->GetBinCenter(KAON_PLUS_BIN));
                }
                particles[j].Decay2Body(particles[N_PARTICLES_PER_ITERATION + 2 * nDecayedParticles], particles[N_PARTICLES_PER_ITERATION + 2 * nDecayedParticles + 1]);
                nDecayedParticles++;
            }

            // Fill generation histograms
            azimutAngleH->Fill(phi);
            polarAngleH->Fill(theta);
            momentumH->Fill(P);
            transverseMomentumH->Fill(P * sin(theta));
            particleEnergyH->Fill(particles[j].TotEnergy());
        }

        // Compute invariant masses and fill histograms
        for (int j = 0; j < N_PARTICLES_PER_ITERATION + 2 * nDecayedParticles; j++) {
            
            // Alias particles[j] as p1
            const ParticleType *p1 = particles[j].GetParticleType();

            // Ignore K* particles for invariant mass histograms (charge == 0)
            if (p1->GetCharge() != 0) {
                
                // Iterate over all previous particles in the array
                for (int k = 0; k < j; k++) {
                    
                    // Alias particles[k] as p2
                    const ParticleType *p2 = particles[k].GetParticleType();
                    
                    // Ignore K* particles
                    if (p2->GetCharge() != 0) {

                        // Compute invariant mass and fill histogram
                        invMass = particles[j].InvMass(&particles[k]);
                        invMassH->Fill(invMass);
                        
                        // Fill discordant or concordant charge histogram
                        if (p1->GetCharge() != p2->GetCharge())
                            discordantInvMassH->Fill(invMass);
                        else
                            concordantInvMassH->Fill(invMass);

                        // Fill discordant or concordant pion-kaon histogram
                        if (     (p1->GetName() == "π+" && p2->GetName() == "K-") ||
                                 (p1->GetName() == "π-" && p2->GetName() == "K+") ||
                                 (p1->GetName() == "K+" && p2->GetName() == "π-") ||
                                 (p1->GetName() == "K-" && p2->GetName() == "π+"))
                            discordantPionKaonInvMassH->Fill(invMass);
                        else if ((p1->GetName() == "π+" && p2->GetName() == "K+") ||
                                 (p1->GetName() == "π-" && p2->GetName() == "K-") ||
                                 (p1->GetName() == "K+" && p2->GetName() == "π+") ||
                                 (p1->GetName() == "K-" && p2->GetName() == "π-"))
                            concordantPionKaonInvMassH->Fill(invMass);

                    }
                }
            }

            // Fill histogram with invariant masses of K* daughters
            if (j >= N_PARTICLES_PER_ITERATION && j % 2 == 0)
                daughtersInvMassH->Fill(particles[j].InvMass(&particles[j + 1]));
        }
    }

    file->Write();
    file->Close();
}