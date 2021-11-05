#include "Particle.h"

#include "ParticleType.h"
#include "ResonanceType.h"
#include <string>
#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace std;

int Particle::fNParticleType = 0;
const ParticleType *Particle::fParticleType[Particle::fMaxNumParticleType];

int Particle::FindParticle(string particleName) {
    for (int i = 0; i < fNParticleType; i++)
        if (fParticleType[i]->GetName() == particleName)
            return i;
    return -1;
}

Particle::Particle() {
    fIndex = -1;
}

Particle::Particle(string particleName, double Px = 0, double Py = 0, double Pz = 0) {
    fIndex = FindParticle(particleName);
    if (fIndex < 0)
        std::cout << "Particle " << particleName << " not found" << std::endl;
    fPx = Px;
    fPy = Py;
    fPz = Pz;
}

int Particle::GetIndex() const {
    return fIndex;
}

const ParticleType *Particle::GetParticleType() const {
    return fParticleType[fIndex];
}

void Particle::SetIndex(int index) {
    if (index >= fNParticleType) {
        std::cout << "Cannot set index " << index << ": out of bounds error" << std::endl;
        return;
    }
    fIndex = index;
}

void Particle::SetIndex(string particleName) {
    int index = FindParticle(particleName);
    if (index < 0) {
        std::cout << "Cannot set index for searched particle " << particleName << ": particle not found" << std::endl;
        return;
    }
    fIndex = index;
}

void Particle::AddParticleType(const string particleName, const double mass, const int charge, const double width) {
    if (fNParticleType >= fMaxNumParticleType) {
        std::cout << "Maximum number of particle types reached: cannot add new particle type" << std::endl;
        return;
    }

    const int index = FindParticle(particleName);
    if (index >= 0) {
        std::cout << "Particle type " << particleName << " already exists: cannot add duplicate" << std::endl;
        return;
    }

    fParticleType[fNParticleType] = width == 0 ? new ParticleType(particleName, mass, charge) : new ResonanceType(particleName, mass, charge, width);
    fNParticleType++;
}

void Particle::PrintParticleTypes() {
    for (int i = 0; i < fNParticleType; i++)
        fParticleType[i]->Print();
}

void Particle::Print() const {
    std::cout << fParticleType[fIndex]->GetName() << " [index = " << fIndex << "]" << std::endl <<
                 "\tP = (" << fPx << ", " << fPy << ", " << fPz << ")" << std::endl;
}

double Particle::GetPx() const {
    return fPx;
}

double Particle::GetPy() const {
    return fPy;
}

double Particle::GetPz() const {
    return fPz;
}

double Particle::GetMass() const {
    return fParticleType[fIndex]->GetMass();
}

double Particle::TotEnergy() const {
    return sqrt(pow(GetMass(), 2) + pow(fPx, 2) + pow(fPy, 2) + pow(fPz, 2));
}

double Particle::InvMass(Particle *p) const {
    return sqrt(pow(TotEnergy() + p->TotEnergy(), 2) - (pow(fPx + p->GetPx(), 2) + pow(fPy + p->GetPy(), 2) + pow(fPz + p->GetPz(), 2)));
}

void Particle::SetP(double Px, double Py, double Pz) {
    fPx = Px;
    fPy = Py;
    fPz = Pz;
}

int Particle::Decay2Body(Particle &dau1, Particle &dau2) const {
    if (GetMass() == 0.0) {
        printf("Decayment cannot be preformed if mass is zero\n");
        return 1;
    }

    double massMot = GetMass();
    double massDau1 = dau1.GetMass();
    double massDau2 = dau2.GetMass();

    // add width effect
    if (fIndex > -1) {
        // gaussian random numbers
        float x1, x2, w, y1, y2;

        double invnum = 1. / RAND_MAX;
        do {
            x1 = 2.0 * rand() * invnum - 1.0;
            x2 = 2.0 * rand() * invnum - 1.0;
            w = x1 * x1 + x2 * x2;
        } while (w >= 1.0);

        w = sqrt((-2.0 * log(w)) / w);
        y1 = x1 * w;
        y2 = x2 * w;

        massMot += fParticleType[fIndex]->GetWidth() * y1;
    }

    if (massMot < massDau1 + massDau2) {
        printf("Decayment cannot be preformed because mass is too low in this channel\n");
        return 2;
    }

    double pout = sqrt((massMot * massMot - (massDau1 + massDau2) * (massDau1 + massDau2)) * (massMot * massMot - (massDau1 - massDau2) * (massDau1 - massDau2))) / massMot * 0.5;

    double norm = 2 * M_PI / RAND_MAX;

    double phi = rand() * norm;
    double theta = rand() * norm * 0.5 - M_PI / 2.;
    dau1.SetP(pout * sin(theta) * cos(phi), pout * sin(theta) * sin(phi), pout * cos(theta));
    dau2.SetP(-pout * sin(theta) * cos(phi), -pout * sin(theta) * sin(phi), -pout * cos(theta));

    double energy = sqrt(fPx * fPx + fPy * fPy + fPz * fPz + massMot * massMot);

    double bx = fPx / energy;
    double by = fPy / energy;
    double bz = fPz / energy;

    dau1.Boost(bx, by, bz);
    dau2.Boost(bx, by, bz);

    return 0;
}

void Particle::Boost(double bx, double by, double bz) {
    double energy = TotEnergy();

    // Boost this Lorentz vector
    double b2 = bx * bx + by * by + bz * bz;
    double gamma = 1.0 / sqrt(1.0 - b2);
    double bp = bx * fPx + by * fPy + bz * fPz;
    double gamma2 = b2 > 0 ? (gamma - 1.0) / b2 : 0.0;

    fPx += gamma2 * bp * bx + gamma * bx * energy;
    fPy += gamma2 * bp * by + gamma * by * energy;
    fPz += gamma2 * bp * bz + gamma * bz * energy;
}