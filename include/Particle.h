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
    /* BEGIN Constructor & Destructor */
    Particle(const double charge, const double mass0, const TVector3 &momentum, const Detector &detector);
    ~Particle();
    /* END Constructor & Destructor */

    /* BEGIN Getters */
    // returns the charge of the particle in e
    double getCharge() const;

    // returns the rest mass of the particle in GeV/c^2
    double getMass0() const;

    // returns the momentum of the particle in GeV/c
    TVector3 getMomentum() const;

    // returns the position of the particle in cm
    TVector3 getPosition() const;
    /* END Getters */
};

#endif /* PARTICLE_H */
