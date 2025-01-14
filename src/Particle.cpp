#include "Particle.h"

#include <TMath.h>

Particle::Particle(const double charge, const double mass0, const TVector3 &momentum, const TVector3 &position)
    : charge(charge), mass0(mass0), momentum(momentum), position(position)
{
    energy = TMath::Sqrt(momentum.Mag2() + mass0 * mass0);
    gamma = energy / mass0;
    mass = gamma * mass0;
    velocity = momentum * (1. / energy);
    beta = velocity.Mag();
}

Particle::Particle(const double charge, const double mass0, double beta, const TVector3 &position, const TVector3 &direction)
    : charge(charge), mass0(mass0), position(position)
{
    if (beta < 0)
    {
        printf("Beta cannot be negative, setting beta to 0!\n");
        beta = 0;
    }

    if (beta >= 1)
    {
        printf("Beta must be less than 1, setting beta to 0.9999999999999999!\n");
        beta = 0.9999999999999999;
    }

    this->beta = beta;
    energy = mass0 / TMath::Sqrt(1 - beta * beta);
    gamma = energy / mass0;
    mass = gamma * mass0;
    velocity = beta * direction;
    momentum = mass * velocity;
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

double Particle::getMass() const
{
    return mass;
}

double Particle::getEnergy() const
{
    return energy;
}

double Particle::getGamma() const
{
    return gamma;
}

double Particle::getBeta() const
{
    return beta;
}

TVector3 Particle::getPosition() const
{
    return position;
}

TVector3 Particle::getMomentum() const
{
    return momentum;
}

TVector3 Particle::getVelocity() const
{
    return velocity;
}

bool Particle::setBeta(const double beta)
{
    if (beta < 0 || beta >= 1)
    {
        printf("Beta must be between 0 and 1!\n");
        return false;
    }

    this->beta = beta;
    energy = mass0 / TMath::Sqrt(1 - beta * beta);
    gamma = energy / mass0;
    mass = gamma * mass0;
    velocity = beta * velocity.Unit();
    momentum = mass * velocity;

    return true;
}

bool Particle::setMomentum(const TVector3 &momentum)
{
    this->momentum = momentum;
    energy = TMath::Sqrt(momentum.Mag2() + mass0 * mass0);
    gamma = energy / mass0;
    mass = gamma * mass0;
    velocity = momentum * (1. / energy);
    beta = velocity.Mag();

    return true;
}

bool Particle::setBetaGamma(const double betaGamma)
{
    if (betaGamma < 0)
    {
        printf("Beta * gamma cannot be negative!\n");
        return false;
    }

    beta = betaGamma / TMath::Sqrt(1 + betaGamma * betaGamma);
    gamma = betaGamma / beta;
    energy = gamma * mass0;
    mass = gamma * mass0;
    velocity = beta * velocity.Unit();
    momentum = mass * velocity;

    return true;
}

double Particle::Wmax() const
{
    constexpr double electronMass = 0.51099895000; // Rest mass of the electron in MeV/c^2

    const double betaGamma = beta * gamma;                     // Beta * gamma of the particle
    const double electronMassMassRatio = electronMass / mass0; // Rest mass of the electron / rest mass of the particle

    return 2 * electronMass * betaGamma * betaGamma / (1 + 2 * gamma * electronMassMassRatio + electronMassMassRatio * electronMassMassRatio);
}
