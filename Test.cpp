#include "Parameters.h"
#include "Particle.h"
#include "ParticleType.h"
#include "ResonanceType.h"

#include <TFile.h>
#include <iostream>

using namespace std;

void Test() {
    TFile *sourceFile = new TFile("histograms.root", "READ");
    TFile *destinationFile = new TFile("analysis.root", "RECREATE");

    cout << objects->GetSize() << endl;
    for (TObject *&&obj : *objects) {
        cout << "ciao" << endl;
        destinationFile->WriteTObject(obj);
    }
    sourceFile->Close();
    delete sourceFile;
    destinationFile->Close();
    delete destinationFile;
}