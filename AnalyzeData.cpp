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

    // Set styles for histograms and fit curves
    gStyle->SetOptStat("e");
    gStyle->SetOptFit(1111);

    particleTypesH->SetFillColor(kBlue);
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

    azimutAngleH->GetYaxis()->SetTitle("Probability Density");
    polarAngleH->GetYaxis()->SetTitle("Probability Density");
    momentumH->GetYaxis()->SetTitle("Probability Density");
    discordantInvMassH->GetYaxis()->SetTitle("Occurrences");
    concordantInvMassH->GetYaxis()->SetTitle("Occurrences");
    discordantPionKaonInvMassH->GetYaxis()->SetTitle("Occurrences");
    concordantPionKaonInvMassH->GetYaxis()->SetTitle("Occurrences");
    daughtersInvMassH->GetYaxis()->SetTitle("Occurrences");

    // Output data for tables
    cout << "Particle type occurrences:" << endl;
    for (int i = 1; i <= N_PARTICLE_TYPES; i++)
        cout << LABELS[i - 1] << " " << particleTypesH->GetBinContent(i) << " +/- " << particleTypesH->GetBinError(i) << endl;

    // Fit angles distribution and check uniformity
    azimutAngleH->Scale(1.0 / azimutAngleH->Integral("width"));
    azimutAngleH->Fit("pol0", "Q");
    polarAngleH->Scale(1.0 / polarAngleH->Integral("width"));
    polarAngleH->Fit("pol0", "Q");

    // Fit momentum with exponential distribution
    TF1 *exponential = new TF1("exponential", "[0]*exp(-[0]*x)", 0.0, TMath::Infinity());
    exponential->SetParameter(0, 1);
    momentumH->Scale(1.0 / momentumH->Integral("width"));
    momentumH->Fit("exponential", "QW");

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

    // Save histograms in files for the final report
    TCanvas *canvas;
    
    canvas = new TCanvas();
    canvas->Divide(2, 2);
    canvas->cd(1);
    particleTypesH->Draw();
    canvas->cd(2);
    azimutAngleH->Draw();
    canvas->cd(3);
    polarAngleH->Draw();
    canvas->cd(4);
    momentumH->Draw();
    canvas->SaveAs("./histograms/typesAndMomentum.tikz.tex");
    delete canvas;

    canvas = new TCanvas();
    canvas->Divide(1, 3);
    canvas->cd(1);
    daughtersInvMassH->Draw();
    canvas->cd(2);
    discordantMinusConcordantH->Draw();
    canvas->cd(3);
    pionKaonDiscordantMinusConcordantH->Draw();
    canvas->SaveAs("./histograms/invMass.tikz.tex");
    delete canvas;

    // Close root file
    file->Close();
    
}