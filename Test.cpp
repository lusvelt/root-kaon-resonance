#include "Parameters.h"
#include "Particle.h"
#include "ParticleType.h"
#include "ResonanceType.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <TH1D.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TROOT.h>

void Test() {
    ofstream output;
    cout.rdbuf(output.rdbuf());
    output.open("output.txt");

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

    // Variables definitions
    Particle particles[N_PARTICLE_TYPES + MAX_PRODUCTS];
    double phi, theta, P, rndm;
    double Px, Py, Pz;
    double invMass;

    int nDecayedParticles; // Counter of decayed particles

    // Iterations
    for (int i = 0; i < 1; i++) {
        nDecayedParticles = 0;
        for (int j = 0; j < N_PARTICLES_PER_ITERATION; j++) {
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
            if (rndm < PION_PLUS_CUMULATIVE)
                particles[j].SetIndex("π+");
            else if (rndm < PION_MINUS_CUMULATIVE)
                particles[j].SetIndex("π-");
            else if (rndm < KAON_PLUS_CUMULATIVE)
                particles[j].SetIndex("K+");
            else if (rndm < KAON_MINUS_CUMULATIVE)
                particles[j].SetIndex("K-");
            else if (rndm < PROTON_PLUS_CUMULATIVE)
                particles[j].SetIndex("p+");
            else if (rndm < PROTON_MINUS_CUMULATIVE)
                particles[j].SetIndex("p-");
            else {
                particles[j].SetIndex("K*");

                // Decayment of K*
                rndm = gRandom->Rndm();
                if (rndm < 0.5) {
                    particles[N_PARTICLES_PER_ITERATION + 2 * nDecayedParticles].SetIndex("π+");
                    particles[N_PARTICLES_PER_ITERATION + 2 * nDecayedParticles + 1].SetIndex("K-");
                } else {
                    particles[N_PARTICLES_PER_ITERATION + 2 * nDecayedParticles].SetIndex("π-");
                    particles[N_PARTICLES_PER_ITERATION + 2 * nDecayedParticles + 1].SetIndex("K+");
                }
                particles[j].Decay2Body(particles[N_PARTICLES_PER_ITERATION + 2 * nDecayedParticles], particles[N_PARTICLES_PER_ITERATION + 2 * nDecayedParticles + 1]);
                nDecayedParticles++;
            }
        }
        
        for (int j = 0; j < N_PARTICLES_PER_ITERATION + 2 * nDecayedParticles; j++) {
            const ParticleType *p1 = particles[j].GetParticleType();
            cout << p1->GetName() << endl;
            // Ignore K* particles for invariant mass histograms (charge == 0)
            if (p1->GetCharge() != 0) {
                for (int k = 0; k < j; k++) {
                    const ParticleType *p2 = particles[k].GetParticleType();
                    cout << "\t" << p2->GetName() << " ";
                    if (p2->GetCharge() != 0) {
                        invMass = particles[j].InvMass(&particles[k]);

                        if (p1->GetCharge() != p2->GetCharge()) {
                            discordantInvMassH->Fill(invMass);
                            cout << "DISCORDANT ";
                        } else {
                            concordantInvMassH->Fill(invMass);
                            cout << "CONC ";
                        }

                        if ((p1->GetName() == "π+" && p2->GetName() == "K-") ||
                            (p1->GetName() == "π-" && p2->GetName() == "K+")) {
                            discordantPionKaonInvMassH->Fill(invMass);
                            cout << "PION/KAON ";
                        } else if ((p1->GetName() == "π+" && p2->GetName() == "K+") ||
                                 (p1->GetName() == "π-" && p2->GetName() == "K-")) {
                            concordantPionKaonInvMassH->Fill(invMass);
                            cout << "PION/KAON ";
                        }
                    }
                    cout << endl;
                }
            }
        }


    }
    output.close();
}