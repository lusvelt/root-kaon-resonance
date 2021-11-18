#include "Parameters.h"

#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>
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
    TH1D *azimutAngleH = (TH1D*) file->Get("azimutAngleH");
    TH1D *polarAngleH = (TH1D*) file->Get("polarAngleH");
    TH1D *momentumH = (TH1D*) file->Get("momentumH");
    TH1D *discordantInvMassH = (TH1D*) file->Get("discordantInvMassH");
    TH1D *concordantInvMassH = (TH1D*) file->Get("concordantInvMassH");
    TH1D *discordantPionKaonInvMassH = (TH1D*) file->Get("discordantPionKaonInvMassH");
    TH1D *concordantPionKaonInvMassH = (TH1D*) file->Get("concordantPionKaonInvMassH");
    TH1D *daughtersInvMassH = (TH1D*) file->Get("daughtersInvMassH");

    // Set styles for histograms
    gStyle->SetOptStat("e");
    gStyle->SetOptFit(1111);
    particleTypesH->SetFillColor(kBlue);

    // Set axes labels
    particleTypesH->GetXaxis()->SetBinLabel(PION_PLUS_BIN, "#pi+");
    particleTypesH->GetXaxis()->SetBinLabel(PION_MINUS_BIN, "#pi-");
    particleTypesH->GetXaxis()->SetBinLabel(KAON_PLUS_BIN, "K+");
    particleTypesH->GetXaxis()->SetBinLabel(KAON_MINUS_BIN, "K-");
    particleTypesH->GetXaxis()->SetBinLabel(PROTON_PLUS_BIN, "p+");
    particleTypesH->GetXaxis()->SetBinLabel(PROTON_MINUS_BIN, "p+");
    particleTypesH->GetXaxis()->SetBinLabel(KAON_STAR_BIN, "K*");
    particleTypesH->GetYaxis()->SetTitle("Occurrences");

    azimutAngleH->GetXaxis()->SetTitle("Angle (rad)");
    polarAngleH->GetXaxis()->SetTitle("Angle (rad)");
    momentumH->GetXaxis()->SetTitle("Momentum (GeV/c)");
    discordantInvMassH->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
    concordantInvMassH->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
    discordantPionKaonInvMassH->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
    concordantPionKaonInvMassH->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
    daughtersInvMassH->GetXaxis()->SetTitle("Mass (GeV/c^{2})");

    azimutAngleH->GetYaxis()->SetTitle("Frequency Density");
    polarAngleH->GetYaxis()->SetTitle("Frequency Density");
    momentumH->GetYaxis()->SetTitle("Frequency Density");
    discordantInvMassH->GetYaxis()->SetTitle("Occurrences");
    concordantInvMassH->GetYaxis()->SetTitle("Occurrences");
    discordantPionKaonInvMassH->GetYaxis()->SetTitle("Occurrences");
    concordantPionKaonInvMassH->GetYaxis()->SetTitle("Occurrences");
    daughtersInvMassH->GetYaxis()->SetTitle("Frequency Density");

    // Output partile type occurrences for the report table
    cout << "Particle type occurrences:" << endl;
    for (int i = 1; i <= N_PARTICLE_TYPES; i++)
        cout << "\t" << LABELS[i - 1] << " " << particleTypesH->GetBinContent(i) << " +/- " << particleTypesH->GetBinError(i) << endl;

    // Convert occurrencies into frequency density
    azimutAngleH->Scale(1.0 / azimutAngleH->Integral("width"));
    polarAngleH->Scale(1.0 / polarAngleH->Integral("width"));
    momentumH->Scale(1.0 / momentumH->Integral("width"));
    daughtersInvMassH->Scale(1.0 / daughtersInvMassH->Integral("width"));
    
    // Fit histograms with the right distribution
    azimutAngleH->Fit("pol0", "Q");         // Uniform
    polarAngleH->Fit("pol0", "Q");          // Uniform
    momentumH->Fit("expo", "Q");            // Exponential
    daughtersInvMassH->Fit("gaus", "Q");    // Normal

    // Output mean for momentum fit distribution
    TF1 *momentumFit = momentumH->GetFunction("expo");
    double momentumFitMean = -1.0 / momentumFit->GetParameter(1);
    double momentumFitMeanError = abs(momentumFitMean * momentumFit->GetParError(1) / momentumFit->GetParameter(1));
    cout << "Momentum Exponential Fit:" << endl;
    cout << "\tMean: " << momentumFitMean << " +/- " << momentumFitMeanError << endl;

    // Check consistency of invariant mass histograms
    TH1D *pionKaonDiscordantMinusConcordantH = (TH1D*) discordantPionKaonInvMassH->Clone("pionKaonDiscordantMinusConcordantH");
    pionKaonDiscordantMinusConcordantH->Add(concordantPionKaonInvMassH, -1.0);
    pionKaonDiscordantMinusConcordantH->Fit("gaus", "Q");
    pionKaonDiscordantMinusConcordantH->SetTitle("Discordant-Concordant Pion/Kaon Invariant Mass Difference");

    TH1D *discordantMinusConcordantH = (TH1D*) discordantInvMassH->Clone("discordantMinusConcordantH");
    discordantMinusConcordantH->Add(concordantInvMassH, -1.0);
    discordantMinusConcordantH->Fit("gaus", "Q");
    discordantMinusConcordantH->SetTitle("Discordant-Concordant Invariant Mass Difference");

    // Save histograms in files for the final report
    TCanvas *c1 = new TCanvas();
    c1->Divide(2, 2);
    c1->cd(1);
    particleTypesH->Draw();
    c1->cd(2);
    azimutAngleH->Draw();
    c1->cd(3);
    polarAngleH->Draw();
    c1->cd(4);
    momentumH->Draw();
    c1->SaveAs("./histograms/typesAndMomentum.tikz.tex");

    gStyle->SetOptStat(0);
    TCanvas *c2 = new TCanvas("c2", "Invariant Mass", 600, 800);
    c2->Divide(1, 3);
    c2->cd(1);
    daughtersInvMassH->Draw();
    c2->cd(2);
    discordantMinusConcordantH->Draw();
    c2->cd(3);
    pionKaonDiscordantMinusConcordantH->Draw();
    c2->SaveAs("./histograms/invMass.tikz.tex");

    // Close root file
    file->Close();
}