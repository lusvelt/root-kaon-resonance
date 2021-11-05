#include "ParticleType.h"

#ifndef PARTICLE_H
#define PARTICLE_H

using namespace std;

class Particle {
    public:
        Particle();
        Particle(string particleName, double Px, double Py, double Pz);
        int GetIndex() const;
        const ParticleType *GetParticleType() const;
        void SetIndex(int index);
        void SetIndex(string particleName);
        void Print() const;
        double GetPx() const;
        double GetPy() const;
        double GetPz() const;
        double GetMass() const;
        double TotEnergy() const;
        double InvMass(Particle *p) const;
        void SetP(double Px, double Py, double Pz);
        int Decay2Body(Particle &dau1, Particle &dau2) const;

        static void AddParticleType(string particleName, const double mass, const int charge, const double width = 0);
        static void PrintParticleTypes();

    private:
        int fIndex;
        double fPx, fPy, fPz;

        static const int fMaxNumParticleType = 10;
        static const ParticleType *fParticleType[];
        static int fNParticleType;

        static int FindParticle(string particleName);

        void Boost(double bx, double by, double bz);
};

#endif