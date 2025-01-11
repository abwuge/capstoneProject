// Coordinate system definition
/******************************\
*           ·------- x         *
*          /|                  *
*         / |                  *
*        y  |                  *
*           z                  *
\******************************/

#include <iostream>
#include <TVector3.h>

#include "Particle.h"
#include "calculator.h"

int main(int argc, char *argv[])
{
    // Particle properties
    const double charge = -1;              // Charge of the particle in e
    const double mass0 = 0.105658;         // Rest mass of the particle in GeV/c^2
    const TVector3 startMomentum(0, 0, 1); // Initial momentum of the particle in GeV/c
    const TVector3 startPosition(0, 0, 0); // Initial position of the particle in cm

    Particle startMuon(charge, mass0, startMomentum, startPosition);

    // Environment properties
    const double c = 299792458;   // speed of light in m/s
    const TVector3 B(0, 0.14, 0); // Magnetic field in Tesla

    const std::vector<double> scintillatorCountersLocation = {-61, -60, 60, 61}; // z-coordinates of the scintillator counters in cm (assuming they are infinite planes in the x and y directions)

    std::vector<TVector3> hitPositions = getHitPositions(startMuon, B, scintillatorCountersLocation);

    return 0;
}