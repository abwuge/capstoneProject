#include <iostream>
#include <TVector3.h>

#include "Material.h"
#include "ScintillatorCounters.h"
#include "Particle.h"
#include "Detector.h"

int main(int argc, char *argv[])
{
    /* BEGIN Material properties */
    // The material is polyvinyltoluene [(2-CH3C6H4CHCH2)n]
    // The data below is taken from https://pdg.lbl.gov/2024/AtomicNuclearProperties/HTML/polyvinyltoluene.html (change "2024" in the URL to "current" for the latest data)
    // Although the density of EJ-200 is 1.023 g/cm^3 (from https://eljentechnology.com/products/plastic-scintillators/ej-200-ej-204-ej-208-ej-212), the density of polyvinyltoluene is used
    constexpr double density = 1.032;                     // Density of polyvinyltoluene in g/cm^3
    constexpr double Z = 1 * 10.00 + 6 * 9.03;            // Atomic number of polyvinyltoluene (H: 10, C: 9.03)
    constexpr double A = 1.0080 * 10.00 + 12.0107 * 9.03; // Atomic mass of polyvinyltoluene (H: 10, C: 9.03)
    constexpr double zToARatio = 0.54141;                 // <Z/A> of polyvinyltoluene in g/mol
    constexpr double I = 64.7;                            // Mean excitation energy of polyvinyltoluene in eV
    constexpr double hbarOmegaP = 21.54;                  // Plasma energy of polyvinyltoluene in eV

    const Material polyvinyltoluene(density, Z, A, zToARatio, I, hbarOmegaP); // Polyvinyltoluene material
    /* END Material properties */

    /* BEGIN Detector properties */
    const std::vector<double> EJ_200Location = {-65, -63, 63, 65}; // z-coordinate of EJ-200 scintillator counters in AMS-02 in cm (assuming they are infinite planes in the x and y directions)
    const std::vector<bool> EJ_200Direction = {1, 0, 0, 1};        // Direction of EJ-200 scintillator counters in AMS-02 (1 for x, 0 for y)
    constexpr double thickness = 1;                                // Thickness of EJ-200 scintillator counters in cm

    std::vector<ScintillatorCounters> EJ_200; // EJ-200 scintillator counters in AMS-02
    EJ_200.reserve(4);
    for (int i = 0; i < 4; i++)
    {
        EJ_200.push_back(ScintillatorCounters(EJ_200Location[i], EJ_200Direction[i], thickness, polyvinyltoluene));
    }
    const TVector3 B(0, 0.14, 0); // Magnetic field in AMS-02 in Tesla

    const Detector TOF(EJ_200, B); // Time-of-flight detector in AMS-02
    // /* END Detector properties */

    // /* BEGIN Particle properties */
    // // The initial position of the particle is fixed at the center of the first scintillator counter
    // constexpr double charge = -1;                      // Charge of muon in e
    // constexpr double mass0 = 105.658;                  // Rest mass of muon in MeV/c^2
    // const TVector3 startMomentum(0, 0, 1e3);           // Initial momentum of muon in MeV/c
    // const TVector3 startPosition(0, 0, TOF.getMinZ()); // Initial position of muon in cm

    // const Particle startMuon(charge, mass0, startMomentum, startPosition); // Initial muon
    // /* END Particle properties */

    /* BEGIN Particle properties */
    // The initial position of the particle is fixed at the center of the first scintillator counter
    constexpr double charge = 3;                       // Charge of Li6 in e
    constexpr double u = 931.49410372;                 // Atomic mass unit in MeV/c^2 (from https://en.wikipedia.org/wiki/Dalton_(unit))
    constexpr double mass0 = 6.0151228874 * u;         // Rest mass of Li6 in MeV/c^2 (from https://en.wikipedia.org/wiki/Isotopes_of_lithium)
    const TVector3 startMomentum(0, 0, 1e4);           // Initial momentum of Li6 in MeV/c
    const TVector3 startPosition(0, 0, TOF.getMinZ()); // Initial position of Li6 in cm

    const Particle startLi6(charge, mass0, startMomentum, startPosition); // Initial Li6
    /* END Particle properties */

    std::cout << "Cyclotron radius: " << TOF.particleCyclotronRadius(startLi6) << " cm" << std::endl;
    std::cout << "beta: " << startLi6.getBeta() << std::endl;
    std::cout << "Stopping power: " << polyvinyltoluene.meanRateOfEnergyLoss(startLi6) << " MeV/cm" << std::endl;

    return 0;
}