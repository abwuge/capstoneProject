#include "Calculator.h"

double Calculator::getCyclotronRadius(const Particle &particle, const Detector &detector)
{
    double absCharge = std::abs(particle.getCharge());
    const TVector3 &momentum = particle.getMomentum();
    const TVector3 &B = detector.getB();

    return momentum.Perp(B) / (absCharge * B.Mag());
}

TVector3 Calculator::getCyclotronDirection(const Particle &particle, const Detector &detector)
{
    double charge = particle.getCharge();
    double chargeSign = charge / std::abs(charge);
    const TVector3 &momentumUnit = particle.getMomentum().Unit();
    const TVector3 &BUint = detector.getB().Unit();

    return momentumUnit.Cross(BUint) * chargeSign;
}

std::vector<TVector3> Calculator::getHitPositions(const Particle &particle, const Detector &detector)
{
    std::vector<TVector3> hitPositions;

    // Particle properties
    const double charge = particle.getCharge();
    const double mass0 = particle.getMass0();
    const TVector3 momentum = particle.getMomentum();
    const TVector3 position = particle.getPosition();

    return hitPositions;
}
