#include "Particle.h"
#include <TMath.h>

Particle::Particle(const double charge, const double mass0, const TVector3 &momentum, const TVector3 &position)
    : charge(charge), mass0(mass0), momentum(momentum), position(position)
{
    energy = TMath::Sqrt(momentum.Mag2() + mass0 * mass0);
    lorentzGamma = energy / mass0;
    mass = mass0 * lorentzGamma;
    velocity = momentum * (1. / energy);
}

Particle::~Particle() {}

double Particle::getCharge() const
{
    return charge;
}

double Particle::getMass0() const
{
    return mass0;
}

TVector3 Particle::getMomentum() const
{
    return momentum;
}

TVector3 Particle::getPosition() const
{
    return position;
}

double Particle::getEnergy() const
{
    return energy;
}

double Particle::getLorentzGamma() const
{
    return lorentzGamma;
}

double Particle::getMass() const
{
    return mass;
}

TVector3 Particle::getVelocity() const
{
    return velocity;
}