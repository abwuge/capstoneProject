#include <iostream>
#include <TVector3.h>

#include "Particle.h"
#include "ScintillatorCounters.h"
#include "Detector.h"
#include "Calculator.h"

int main(int argc, char *argv[])
{
    // Detector properties
    const std::vector<double> scintillatorCountersLocation = {-61, -60, 60, 61}; // z-coordinate of the scintillator counters in cm (assuming they are infinite planes in the x and y directions)
    const std::vector<bool> scintillatorCountersDirection = {1, 0, 0, 1}; // Direction of the scintillator counters
    std::vector<ScintillatorCounters> scintillatorCounters;
    for (int i = 0; i < 4; i++)
    {
        scintillatorCounters.push_back(ScintillatorCounters(scintillatorCountersLocation[i], scintillatorCountersDirection[i]));
    }
    const TVector3 B(0, 0.14, 0); // Magnetic field in Tesla

    const Detector TOF(scintillatorCounters, B);

    // Particle properties
    // The initial position of the particle is fixed at the center of the first scintillator counter
    const double charge = -1;              // Charge of the particle in e
    const double mass0 = 0.105658;         // Rest mass of the particle in GeV/c^2
    const TVector3 startMomentum(0, 0, 1); // Initial momentum of the particle in GeV/c

    const Particle startMuon(charge, mass0, startMomentum, TOF);

    std::vector<TVector3> hitPositions = Calculator::getHitPositions(startMuon, B, scintillatorCountersLocation);

    return 0;
}