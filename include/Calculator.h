#ifndef CALCULATOR_H
#define CALCULATOR_H

#include <vector>
#include <TVector3.h>

#include "Particle.h"

namespace Calculator
{
    // Calculate the cyclotron radius of the particle in the magnetic field
    double getCyclotronRadius(const Particle &particle, const Detector &detector);

    // Calculate the derection of the cyclotron radius
    TVector3 getCyclotronDirection(const Particle &particle, const Detector &detector);

    // Calculate the hit positions of the particle in the scintillator counters
    std::vector<TVector3> getHitPositions(const Particle &particle, const Detector &detector);
}

#endif /* CALCULATOR_H */
