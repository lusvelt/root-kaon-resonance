#include "ParticleType.h"

#ifndef RESONANCE_TYPE_H
#define RESONANCE_TYPE_H

using namespace std;

class ResonanceType : public ParticleType {
    public:
        ResonanceType(const string name, const double mass, const int charge, const double width) :
            ParticleType(name, mass, charge), fWidth(width) {}
        double GetWidth() const;
        void Print() const;
        
    private:
        const double fWidth;
};

#endif