#include "Calculator.h"
#include <TMath.h>

double c = TMath::C();

std::vector<TVector3> Calculator::getHitPositions(const Particle &particle, const TVector3 &B, const std::vector<double> &scintillatorCountersLocation)
{
    std::vector<TVector3> hitPositions;

    // Particle properties
    const double charge = particle.getCharge();
    const double mass0 = particle.getMass0();
    const TVector3 momentum = particle.getMomentum();
    const TVector3 position = particle.getPosition();


    
    return hitPositions;
}
