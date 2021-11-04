#include "ResonanceType.h"

#include "ParticleType.h"
#include <iostream>

double ResonanceType::GetWidth() const {
    return fWidth;
}

void ResonanceType::Print() const {
    ParticleType::Print();
    std::cout << "\tWidth = " << fWidth << " s^-1" << std::endl;
}