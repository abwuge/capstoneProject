#include "calculator.h"

std::vector<TVector3> getHitPositions(Particle& particle, const TVector3& B, const std::vector<double>& scintillatorCountersLocation)
{
    std::vector<TVector3> hitPositions;

    // Particle properties
    const double charge = particle.getCharge();
    const double mass = particle.getMass();
    const TVector3 momentum = particle.getMomentum();
    const TVector3 position = particle.getPosition();


    
    return hitPositions;
}