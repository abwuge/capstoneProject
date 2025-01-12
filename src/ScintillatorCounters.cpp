#include "ScintillatorCounters.h"

ScintillatorCounters::ScintillatorCounters(const double location, const bool direction, const double thickness, const Material &material)
    : location(location), direction(direction), thickness(thickness), material(material) {}

ScintillatorCounters::~ScintillatorCounters() {}

double ScintillatorCounters::getLocation() const
{
    return location;
}

bool ScintillatorCounters::getDirection() const
{
    return direction;
}

double ScintillatorCounters::getThickness() const
{
    return thickness;
}

Material ScintillatorCounters::getMaterial() const
{
    return material;
}

double ScintillatorCounters::energyLoss(const Particle &particle) const
{
    return material.meanRateOfEnergyLoss(particle) * thickness;
}
