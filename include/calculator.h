#ifndef CALCULATOR_H
#define CALCULATOR_H

#include <iostream>
#include <vector>
#include <TVector3.h>

#include "Particle.h"

std::vector<TVector3> getHitPositions(Particle& particle, const TVector3& B, const std::vector<double>& scintillatorCountersLocation);

#endif /* CALCULATOR_H */
