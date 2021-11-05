#include "Parameters.h"

#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TObject.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <cmath>

using namespace std;

void AnalyzeData() {
    gROOT->SetBatch();

    // Open root file and retrieve histograms
    TFile *file = new TFile("histograms.root", "READ");

    TH1D *particleTypesH = (TH1D*) file->Get("particleTypesH");
    TH1D *finalParticleTypesH = (TH1D*) file->Get("finalParticleTypesH");
    TH1D *azimutAngleH = (TH1D*) file->Get("azimutAngleH");
    TH1D *polarAngleH = (TH1D*) file->Get("polarAngleH");
    TH1D *momentumH = (TH1D*) file->Get("momentumH");
    TH1D *transverseMomentumH = (TH1D*) file->Get("transverseMomentumH");
    TH1D *particleEnergyH = (TH1D*) file->Get("particleEnergyH");
    TH1D *invMassH = (TH1D*) file->Get("invMassH");
    TH1D *discordantInvMassH = (TH1D*) file->Get("discordantInvMassH");
    TH1D *concordantInvMassH = (TH1D*) file->Get("concordantInvMassH");
    TH1D *discordantPionKaonInvMassH = (TH1D*) file->Get("discordantPionKaonInvMassH");
    TH1D *concordantPionKaonInvMassH = (TH1D*) file->Get("concordantPionKaonInvMassH");
    TH1D *daughtersInvMassH = (TH1D*) file->Get("daughtersInvMassH");

    // Set styles for histograms and fit curves
    gStyle->SetOptStat("RM");
    gStyle->SetOptFit(111);
    gROOT->ForceStyle();

    // Check number of entries of generation histograms
    if (particleTypesH->GetEntries() != N_PARTICLES_PER_ITERATION * N_ITERATIONS)
        cout << "Number of entries of Particle Types Histogram is incorrect" << endl;
    if (azimutAngleH->GetEntries() != N_PARTICLES_PER_ITERATION * N_ITERATIONS)
        cout << "Number of entries of Azimut Angle Histogram is incorrect" << endl;
    if (polarAngleH->GetEntries() != N_PARTICLES_PER_ITERATION * N_ITERATIONS)
        cout << "Number of entries of Polar Angle Histogram is incorrect" << endl;
    if (momentumH->GetEntries() != N_PARTICLES_PER_ITERATION * N_ITERATIONS)
        cout << "Number of entries of Momentum Histogram is incorrect" << endl;
    if (transverseMomentumH->GetEntries() != N_PARTICLES_PER_ITERATION * N_ITERATIONS)
        cout << "Number of entries of Transverse Momentum Histogram is incorrect" << endl;
    if (particleEnergyH->GetEntries() != N_PARTICLES_PER_ITERATION * N_ITERATIONS)
        cout << "Number of entries of Particle Energy Histogram is incorrect" << endl;

    // Alias number of particles for each type and respective errors
    const double nPionPlus = finalParticleTypesH->GetBinContent(PION_PLUS_BIN);
    const double nPionMinus = finalParticleTypesH->GetBinContent(PION_MINUS_BIN);
    const double nKaonPlus = finalParticleTypesH->GetBinContent(KAON_PLUS_BIN);
    const double nKaonMinus = finalParticleTypesH->GetBinContent(KAON_MINUS_BIN);
    const double nProtonPlus = finalParticleTypesH->GetBinContent(PROTON_PLUS_BIN);
    const double nProtonMinus = finalParticleTypesH->GetBinContent(PROTON_MINUS_BIN);
    const double nKaonStar = finalParticleTypesH->GetBinContent(KAON_STAR_BIN);

    const double nPionPlusErr = finalParticleTypesH->GetBinError(PION_PLUS_BIN);
    const double nPionMinusErr = finalParticleTypesH->GetBinError(PION_MINUS_BIN);
    const double nKaonPlusErr = finalParticleTypesH->GetBinError(KAON_PLUS_BIN);
    const double nKaonMinusErr = finalParticleTypesH->GetBinError(KAON_MINUS_BIN);
    const double nProtonPlusErr = finalParticleTypesH->GetBinError(PROTON_PLUS_BIN);
    const double nProtonMinusErr = finalParticleTypesH->GetBinError(PROTON_MINUS_BIN);
    const double nKaonStarErr = finalParticleTypesH->GetBinError(KAON_STAR_BIN);

    const double nPionPlusRelativeErr = nPionPlusErr / nPionPlus;
    const double nPionMinusRelativeErr = nPionMinusErr / nPionMinus;
    const double nKaonPlusRelativeErr = nKaonPlusErr / nKaonPlus;
    const double nKaonMinusRelativeErr = nKaonMinusErr / nKaonMinus;
    const double nProtonPlusRelativeErr = nProtonPlusErr / nProtonPlus;
    const double nProtonMinusRelativeErr = nProtonMinusErr / nProtonMinus;
    const double nKaonStarRelativeErr = nKaonStarErr / nKaonStar;

    // Compute derived stats of combined particles and respective errors
    const double nFinalParticlesPerIteration = (N_PARTICLES_PER_ITERATION * N_ITERATIONS + nKaonStar) / N_ITERATIONS;
    const double nFinalParticlesPerIterationErr = nKaonStarErr / N_ITERATIONS;
    const double nFinalParticlesPerIterationRelativeErr = nFinalParticlesPerIterationErr / nFinalParticlesPerIteration;

    const double nPairs = (nFinalParticlesPerIteration * (nFinalParticlesPerIteration - 1) / 2) * N_ITERATIONS;
    const double nPairsErr = 2 * nFinalParticlesPerIterationErr / nFinalParticlesPerIteration * nPairs;
    const double nPairsRelativeErr = nPairsErr / nPairs;
    
    const double nPositiveParticlesPerIteration = (nPionPlus + nKaonPlus + nProtonPlus) / N_ITERATIONS;
    const double nPositiveParticlesPerIterationErr = (nPionPlusErr + nKaonPlusErr + nProtonPlusErr) / N_ITERATIONS;
    const double nPositiveParticlesPerIterationRelativeErr = nPositiveParticlesPerIterationErr / nPositiveParticlesPerIteration;
    
    const double nNegativeParticlesPerIteration = (nPionMinus + nKaonMinus + nProtonMinus) / N_ITERATIONS;
    const double nNegativeParticlesPerIterationErr = (nPionMinusErr + nKaonMinusErr + nProtonMinusErr) / N_ITERATIONS;
    const double nNegativeParticlesPerIterationRelativeErr = nNegativeParticlesPerIterationErr / nNegativeParticlesPerIteration;
    
    const double nDiscordantPairs = nPositiveParticlesPerIteration * nNegativeParticlesPerIteration * N_ITERATIONS;
    const double nDiscordantPairsRelativeErr = nPositiveParticlesPerIterationRelativeErr + nNegativeParticlesPerIterationRelativeErr;
    const double nDiscordantPairsErr = nDiscordantPairsRelativeErr * nDiscordantPairs;

    const double nConcordantPairs = ((nPositiveParticlesPerIteration * (nPositiveParticlesPerIteration - 1) / 2) +
                                     (nNegativeParticlesPerIteration * (nNegativeParticlesPerIteration - 1) / 2)) *
                                    N_ITERATIONS;
    const double nConcordantPairsRelativeErr = 2 * (nPositiveParticlesPerIterationRelativeErr + nNegativeParticlesPerIterationRelativeErr);
    const double nConcordantPairsErr = nConcordantPairsRelativeErr * nConcordantPairs;

    const double nPositivePionsPerIteration = nPionPlus / N_ITERATIONS;
    const double nPositivePionsPerIterationErr = nPionPlusErr / N_ITERATIONS;
    const double nPositivePionsPerIterationRelativeErr = nPionPlusRelativeErr;

    const double nNegativePionsPerIteration = nPionMinus / N_ITERATIONS;
    const double nNegativePionsPerIterationErr = nPionMinusErr / N_ITERATIONS;
    const double nNegativePionsPerIterationRelativeErr = nPionMinusRelativeErr;
    
    const double nPositiveKaonsPerIteration = nKaonPlus / N_ITERATIONS;
    const double nPositiveKaonsPerIterationErr = nKaonPlusErr / N_ITERATIONS;
    const double nPositiveKaonsPerIterationRelativeErr = nKaonPlusRelativeErr;

    const double nNegativeKaonsPerIteration = nKaonMinus / N_ITERATIONS;
    const double nNegativeKaonsPerIterationErr = nKaonMinusErr / N_ITERATIONS;
    const double nNegativeKaonsPerIterationRelativeErr = nKaonMinusRelativeErr;

    const double nDiscordantPionKaonPairs = (nPositivePionsPerIteration * nNegativeKaonsPerIteration / 2 +
                                             nNegativePionsPerIteration * nPositiveKaonsPerIteration / 2) *
                                            N_ITERATIONS;
    const double nPositivePionNegativeKaonPairsRelativeErr = nPionPlusRelativeErr + nKaonMinusRelativeErr;
    const double nPositivePionNegativeKaonPairsErr = nPositivePionNegativeKaonPairsRelativeErr * nPositivePionsPerIteration * nNegativeKaonsPerIteration;
    const double nNegativePionPositiveKaonPairsRelativeErr = nPionMinusRelativeErr + nKaonPlusRelativeErr;
    const double nNegativePionPositiveKaonPairsErr = nNegativePionPositiveKaonPairsRelativeErr * nNegativePionsPerIteration * nPositiveKaonsPerIteration;
    const double nDiscordantPionKaonPairsErr = (nPositivePionNegativeKaonPairsErr + nNegativePionPositiveKaonPairsErr) * N_ITERATIONS;

    const double nConcordantPionKaonPairs = (nPositivePionsPerIteration * nPositiveKaonsPerIteration / 2 +
                                             nNegativePionsPerIteration * nNegativeKaonsPerIteration / 2) *
                                            N_ITERATIONS;
    const double nPositivePionPositiveKaonPairsRelativeErr = nPionPlusRelativeErr + nKaonPlusRelativeErr;
    const double nPositivePionPositiveKaonPairsErr = nPositivePionPositiveKaonPairsRelativeErr * nPositivePionsPerIteration * nPositiveKaonsPerIteration;
    const double nNegativePionNegativeKaonPairsRelativeErr = nPionMinusRelativeErr + nKaonMinusRelativeErr;
    const double nNegativePionNegativeKaonPairsErr = nNegativePionNegativeKaonPairsRelativeErr * nNegativePionsPerIteration * nNegativeKaonsPerIteration;
    const double nConcordantPionKaonPairsErr = (nPositivePionPositiveKaonPairsErr + nNegativePionNegativeKaonPairsErr) * N_ITERATIONS;

    const double nTotParticles = nPionPlus + nPionMinus + nKaonPlus + nKaonMinus + nProtonPlus + nProtonMinus + nKaonStar;
    const double nTotParticlesErr = nPionPlusErr + nPionMinusErr + nKaonPlusErr + nKaonMinusErr + nProtonPlusErr + nProtonMinusErr + nKaonStarErr;

    const double nDaughters = nTotParticles - N_PARTICLES_PER_ITERATION * N_ITERATIONS;
    const double nDaughtersErr = nTotParticlesErr;

    const double nDaughterPairs = nDaughters / 2;
    const double nDaughterPairsErr = nDaughtersErr / 2;

    // Check number of entries of invariant mass histograms
    if (abs(invMassH->GetEntries() - nPairs) > ERROR_FACTOR * nPairsErr)
        cout << "Number of entries of Invariant Mass Histogram is incorrect" << endl;

    if (abs(discordantInvMassH->GetEntries() - nDiscordantPairs) > ERROR_FACTOR * nDiscordantPairsErr)
        cout << "Number of entries of discordant invariant mass histogram is incorrect" << endl;

    if (abs(concordantInvMassH->GetEntries() - nConcordantPairs) > ERROR_FACTOR * nConcordantPairsErr)
        cout << "Number of entries of concordant invariant mass histogram is incorrect" << endl;

    if (abs(discordantPionKaonInvMassH->GetEntries() - nDiscordantPionKaonPairs) > ERROR_FACTOR * nDiscordantPionKaonPairsErr)
        cout << "Number of entries of discordant pion/kaon invariant mass histogram is incorrect" << endl;
    
    if (abs(concordantPionKaonInvMassH->GetEntries() - nConcordantPionKaonPairs) > ERROR_FACTOR * nConcordantPionKaonPairsErr)
        cout << "Number of entries of concordant pion/kaon invariant mass histogram is incorrect" << endl;

    if (abs(daughtersInvMassH->GetEntries() - nDaughterPairs) > ERROR_FACTOR * nDaughterPairsErr)
        cout << "Number of entries of daughters invariant mass histogram is incorrect" << endl;

    // Check proportions of particle types
    for (int i = 1; i <= particleTypesH->GetNbinsX(); i++) {
        const double nProb = N_PARTICLES_PER_ITERATION * N_ITERATIONS * PROBABILITIES[i - 1];
        const double nParticles = particleTypesH->GetBinContent(i);
        const double nParticlesErr = particleTypesH->GetBinError(i);
        if (nProb < nParticles - ERROR_FACTOR * nParticlesErr || nProb > nParticles + ERROR_FACTOR * nParticlesErr)
            cout << "Number of " << LABELS[i - 1] << " particles is incorrect" << endl;
    }

    // Fit angles distribution and check uniformity
    azimutAngleH->UseCurrentStyle();
    azimutAngleH->Scale(1.0 / azimutAngleH->Integral(), "width"); // Normalization
    azimutAngleH->Fit("pol0", "Q");

    polarAngleH->Scale(1.0 / polarAngleH->Integral(), "width"); // Normalization
    polarAngleH->Fit("pol0", "Q");

    // Fit momentum with exponential distribution
    momentumH->Fit("expo", "Q");

    // Fit daughter invariant mass with normal distribution
    daughtersInvMassH->Fit("gaus", "Q");

    // Check consistency of invariant mass histograms
    TH1D *pionKaonDiscordantMinusConcordantH = (TH1D *)discordantPionKaonInvMassH->Clone("pionKaonDiscordantMinusConcordantH");
    pionKaonDiscordantMinusConcordantH->Add(concordantPionKaonInvMassH, -1.0);
    pionKaonDiscordantMinusConcordantH->Fit("gaus", "Q");
    pionKaonDiscordantMinusConcordantH->SetTitle("Discordant-Concordant Pion/Kaon Invariant Mass Difference");

    TH1D *discordantMinusConcordantH = (TH1D *)discordantInvMassH->Clone("discordantMinusConcordantH");
    discordantMinusConcordantH->Add(concordantInvMassH, -1.0);
    discordantMinusConcordantH->Fit("gaus", "Q");
    discordantMinusConcordantH->SetTitle("Discordant-Concordant Invariant Mass Difference");

    // Print histograms in files
    TList *histograms = new TList();
    histograms->Add(particleTypesH);
    histograms->Add(azimutAngleH);
    histograms->Add(polarAngleH);
    histograms->Add(momentumH);
    histograms->Add(transverseMomentumH);
    histograms->Add(particleEnergyH);
    histograms->Add(invMassH);
    histograms->Add(discordantInvMassH);
    histograms->Add(concordantInvMassH);
    histograms->Add(discordantPionKaonInvMassH);
    histograms->Add(concordantPionKaonInvMassH);
    histograms->Add(daughtersInvMassH);
    histograms->Add(pionKaonDiscordantMinusConcordantH);
    histograms->Add(discordantMinusConcordantH);

    TCanvas *canvas = new TCanvas();
    for (TObject* &&histogram : *histograms) {
        histogram->Draw();
        string directory = "./histograms/";
        string name = histogram->GetName();
        string extension = ".png";
        string fileName = directory + name + extension;
        canvas->Print(fileName.c_str());
    }

    // Close root file
    file->Close();
}