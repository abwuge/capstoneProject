#include "Particle.h"
#include <TMath.h>

Particle::Particle(const double charge, const double mass0, const TVector3 &momentum, const TVector3 &position)
    : charge(charge), mass0(mass0), momentum(momentum), position(position)
{
    energy = TMath::Sqrt(momentum.Mag2() + mass0 * mass0);
    gamma = energy / mass0;
    mass = mass0 * gamma;
    velocity = momentum * (1. / energy);
    beta = velocity.Mag();
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

double Particle::getGamma() const
{
    return gamma;
}

double Particle::getMass() const
{
    return mass;
}

TVector3 Particle::getVelocity() const
{
    return velocity;
}

double Particle::getBeta() const
{
    return beta;
}

double Particle::Wmax() const
{
    constexpr double electronMass = 0.51099895000; // Rest mass of the electron in MeV/c^2

    const double betaGamma = beta * gamma;                     // Beta * gamma of the particle
    const double electronMassMassRatio = electronMass / mass0; // Rest mass of the electron / rest mass of the particle

    return 2 * electronMass * betaGamma * betaGamma / (1 + 2 * gamma * electronMassMassRatio + electronMassMassRatio * electronMassMassRatio);
}
