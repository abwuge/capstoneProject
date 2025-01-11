#ifndef PARTICLE_H
#define PARTICLE_H

#include <iostream>
#include <vector>
#include <TVector3.h>

#include "Detector.h"

class Particle
{
private:
    double charge;     // Charge of the particle in e
    double mass0;      // Rest mass of the particle in GeV/c^2
    TVector3 momentum; // Momentum of the particle in GeV/c
    TVector3 position; // Position of the particle in cm

public:

    // Constructor & Destructor
    Particle(const double charge, const double mass0, const TVector3 &momentum, const Detector &detector);
    ~Particle();

    // Getters
    double getCharge() const;
    double getMass0() const;
    TVector3 getMomentum() const;
    TVector3 getPosition() const;
};

#endif /* PARTICLE_H */