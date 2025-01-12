#include "Detector.h"

#include <algorithm>

#include <TMath.h>

Detector::Detector(const std::vector<ScintillatorCounters> &scintillatorCounters, const TVector3 &B)
    : scintillatorCounters(scintillatorCounters), B(B)
{
    std::sort(this->scintillatorCounters.begin(), this->scintillatorCounters.end(),
              [](const ScintillatorCounters &a, const ScintillatorCounters &b)
              { return a.getLocation() < b.getLocation(); });
}

Detector::~Detector() {}

double Detector::getMinZ() const
{
    return scintillatorCounters.front().getLocation();
}

std::vector<ScintillatorCounters> Detector::getScintillatorCounters() const
{
    return scintillatorCounters;
}

TVector3 Detector::getB() const
{
    return B;
}

double Detector::particleCyclotronRadius(const Particle &particle) const
{
    double absCharge = std::abs(particle.getCharge());
    const TVector3 &momentum = particle.getMomentum();

    if (B.Mag() < 1e-100)
        return 1e100;

    double radius = momentum.Perp(B) / (absCharge * B.Mag());   // cyclotron radius in MeV / c / (e * T) = 1e6 J / (c * T)
    constexpr double conversionFactor = 1e6 / TMath::C() * 1e2; // conversion factor from MeV / c / (e * T) to cm
    return radius * conversionFactor;
}

TVector3 Detector::particleCyclotronDirection(const Particle &particle) const
{
    double charge = particle.getCharge();
    double chargeSign = charge / std::abs(charge);
    const TVector3 &momentumUnit = particle.getMomentum().Unit();
    const TVector3 &BUint = B.Unit();

    return momentumUnit.Cross(BUint) * chargeSign;
}

std::vector<TVector3> Detector::particleHitPositions(const Particle &particle) const
{
    std::vector<TVector3> hitPositions;

    const double radius = particleCyclotronRadius(particle);
    const TVector3 radiusVector = radius * particleCyclotronDirection(particle);
    const TVector3 center = particle.getPosition() + radiusVector;

    return hitPositions;
}
