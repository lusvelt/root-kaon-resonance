#ifndef PARTICLE_TYPE_H
#define PARTICLE_TYPE_H

class ParticleType {
    public:
        ParticleType(string name, const double mass, const int charge) :
            fName(name), fMass(mass), fCharge(charge) {}
        string GetName() const;
        double GetMass() const;
        int GetCharge() const;
        virtual double GetWidth() const;
        virtual void Print() const;

    protected:
        string fName;
        const double fMass;
        const int fCharge;
};

#endif