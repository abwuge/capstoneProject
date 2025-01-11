#include "Particle.h"

Particle::Particle(double charge, double mass, TVector3 momentum, TVector3 position)
{
    this->charge = charge;
    this->mass = mass;
    this->momentum = momentum;
    this->position = position;
}

Particle::~Particle()
{
}

double Particle::getCharge() const
{
    return this->charge;
}

double Particle::getMass() const
{
    return this->mass;
}

TVector3 Particle::getMomentum() const
{
    return this->momentum;
}

TVector3 Particle::getPosition() const
{
    return this->position;
}