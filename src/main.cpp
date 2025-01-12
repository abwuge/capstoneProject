#include <iostream>
#include <TVector3.h>

#include "Particle.h"
#include "ScintillatorCounters.h"
#include "Detector.h"
#include "Calculator.h"

int main(int argc, char *argv[])
{
    /* BEGIN Detector properties */
    const std::vector<double> scintillatorCountersLocation = {-61, -60, 60, 61}; // z-coordinate of the scintillator counters in cm (assuming they are infinite planes in the x and y directions)
    const std::vector<bool> scintillatorCountersDirection = {1, 0, 0, 1};        // Direction of the scintillator counters

    std::vector<ScintillatorCounters> scintillatorCounters;
    for (int i = 0; i < 4; i++)
    {
        scintillatorCounters.push_back(ScintillatorCounters(scintillatorCountersLocation[i], scintillatorCountersDirection[i]));
    }
    const TVector3 B(0, 0.14, 0); // Magnetic field in Tesla

    const Detector TOF(scintillatorCounters, B); // Time-of-flight detector
    /* END Detector properties */

    /* BEGIN Particle properties */
    // The initial position of the particle is fixed at the center of the first scintillator counter
    const double charge = -1;                          // Charge of the particle in e
    const double mass0 = 0.105658;                     // Rest mass of the particle in GeV/c^2
    const TVector3 startMomentum(0, 0, 1);             // Initial momentum of the particle in GeV/c
    const TVector3 startPosition(0, 0, TOF.getMinZ()); // Initial position of the particle in cm

    const Particle startMuon(charge, mass0, startMomentum, startPosition); // Initial muon
    /* END Particle properties */

    std::cout << "Cyclotron radius: " << Calculator::getCyclotronRadius(startMuon, TOF) << " m" << std::endl;
    std::cout << "Cyclotron direction: " << std::endl
              << "X: " << Calculator::getCyclotronDirection(startMuon, TOF).X() << std::endl
              << "Y: " << Calculator::getCyclotronDirection(startMuon, TOF).Y() << std::endl
              << "Z: " << Calculator::getCyclotronDirection(startMuon, TOF).Z() << std::endl;

    return 0;
}