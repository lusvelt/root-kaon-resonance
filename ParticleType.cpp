#include "ParticleType.h"

#include <iostream>

using namespace std;

string ParticleType::GetName() const {
    return fName;
}

double ParticleType::GetMass() const {
    return fMass;
}

int ParticleType::GetCharge() const {
    return fCharge;
}

double ParticleType::GetWidth() const {
    return 0;
}

void ParticleType::Print() const {
    std::cout << "Particle Type: " << fName << std::endl <<
                 "\tMass = " << fMass << " MeV/c^2" << std::endl <<
                 "\tCharge = " << fCharge << " e" << std::endl;
}