#include "Particle.h"

Particle::Particle(const double charge, const double mass0, const TVector3 &momentum, const Detector &detector)
{
    this->charge = charge;
    this->mass0 = mass0;
    this->momentum = momentum;
    this->position = TVector3(0, 0, detector.getMinZ());
}

Particle::~Particle() {}

double Particle::getCharge() const
{
    return this->charge;
}

double Particle::getMass0() const
{
    return this->mass0;
}

TVector3 Particle::getMomentum() const
{
    return this->momentum;
}

TVector3 Particle::getPosition() const
{
    return this->position;
}
