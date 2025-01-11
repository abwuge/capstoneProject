#ifndef CALCULATOR_H
#define CALCULATOR_H

#include <vector>
#include <TVector3.h>

#include "Particle.h"

namespace Calculator
{
    // Calculate the hit positions of the particle in the scintillator counters
    std::vector<TVector3> getHitPositions(const Particle &particle, const TVector3 &B, const std::vector<double> &scintillatorCountersLocation);

}

#endif /* CALCULATOR_H */
