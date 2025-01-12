#ifndef PARTICLE_H
#define PARTICLE_H

#include <iostream>
#include <vector>
#include <TVector3.h>

#include "Detector.h"

class Particle
{
private:
    double charge;       // Charge of the particle in e
    double mass;         // Mass of the particle in GeV/c^2
    double mass0;        // Rest mass of the particle in GeV/c^2
    double energy;       // Energy of the particle in GeV
    double lorentzGamma; // Lorentz factor of the particle
    TVector3 position;   // Position of the particle in cm
    TVector3 momentum;   // Momentum of the particle in GeV/c
    TVector3 velocity;   // Velocity of the particle in c

public:
    /* BEGIN Constructor & Destructor */
    // uses detector to set the initial position of the particle, others are the properties of the particle
    Particle(const double charge, const double mass0, const TVector3 &momentum, const TVector3 &position);
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

    // returns the energy of the particle in GeV
    double getEnergy() const;

    // returns the Lorentz factor of the particle
    double getLorentzGamma() const;

    // returns the mass of the particle in GeV/c^2
    double getMass() const;

    // returns the velocity of the particle
    TVector3 getVelocity() const;
    /* END Getters */
};

#endif /* PARTICLE_H */
