#include <iostream>
#include <algorithm>
#include <thread>
#include <tuple>
#include <vector>

#include <TROOT.h>
#include <TVector3.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TAxis.h>
#include <TLegend.h>

#include "Config.h"
#include "ThreadPool.h"
#include "Material.h"
#include "ScintillatorCounters.h"
#include "Particle.h"
#include "Detector.h"

int main(int argc, char *argv[])
{
    gROOT->SetStyle("Pub");

#if configEnableMultiThreading
    ROOT::EnableThreadSafety();

    if (argc > 1)
        ThreadPool::getInstance(std::min(std::thread::hardware_concurrency(), (unsigned int)std::stoi(argv[1])));
    else
        ThreadPool::getInstance(std::max(std::thread::hardware_concurrency() - 2, 1u));
    
    printf("[Info] Number of threads: %lu\n", ThreadPool::getInstance().getNumThreads());
#endif

/* BEGIN Material properties */
#if configTodo // IDK some porperties of polyvinyltoluene, so I used polystyrene to test the code
    const Material polyvinyltoluene(materials.at(MaterialName::Polystyrene));
#else
    // The EJ-200 scintillator counter is made of polyvinyltoluene [(2-CH3C6H4CHCH2)n]
    const Material polyvinyltoluene(materials.at(MaterialName::Polyvinyltoluene)); // Polyvinyltoluene material
#endif
    /* END Material properties */

    /* BEGIN Detector properties */
    const std::vector<double> EJ_200Location = {-65, -63, 63, 65}; // z-coordinate of EJ-200 scintillator counters in AMS-02 in cm (assuming they are infinite planes in the x and y directions)
    const std::vector<bool> EJ_200Direction = {1, 0, 0, 1};        // Direction of EJ-200 scintillator counters in AMS-02 (1 for x, 0 for y)
    constexpr double thickness = 1;                                // Thickness of EJ-200 scintillator counters in cm
    constexpr double timeResolution = 0.1;                         // Time resolution of EJ-200 scintillator counters in ns

    std::vector<ScintillatorCounters> EJ_200; // EJ-200 scintillator counters in AMS-02
    EJ_200.reserve(4);
    for (int i = 0; i < 4; i++)
    {
        EJ_200.push_back(ScintillatorCounters(EJ_200Location[i], EJ_200Direction[i], thickness, timeResolution, polyvinyltoluene));
    }
    // const TVector3 B(0, 0.14, 0); // Magnetic field in AMS-02 in Tesla
    const TVector3 B(0, 0, 0); // used for testing

    const Detector TOF(EJ_200, B); // Time-of-flight detector in AMS-02
    Detector TOF2(EJ_200, B);
    TOF2 = TOF;
    /* END Detector properties */

    /* BEGIN Particle properties */
    // The initial position of the particle is fixed at the center of the first scintillator counter
    constexpr double charge = 3;                                                 // Charge of Li6 in e
    constexpr double uIn_kg = 1.66053906660e-27;                                 // Atomic mass unit in kg (from https://www.bipm.org/documents/20126/41483022/SI-Brochure-9.pdf)
    constexpr double u = uIn_kg * TMath::C() * TMath::C() / (1e6 * TMath::Qe()); // Atomic mass unit in MeV/c^2
    constexpr double mass0 = 6.01512289 * u;                                     // Rest mass of Li6 in MeV/c^2 (from https://ciaaw.org/lithium.htm)
    constexpr double startBeta = 0.4;                                            // Initial beta of Li6
    const TVector3 startPosition(0, 0, TOF.getMinZ());                           // Initial position of Li6 in cm

    Particle Li6(charge, mass0, startBeta, startPosition); // Initial Li6
    /* END Particle properties */

    // EJ_200.at(0).plotEnergyLoss(Li6, 0.1, 1000, 1000, true, "test/energyLoss.png");
    // EJ_200.at(0).plotEnergyLossFluctuation(Li6, 100000, true, "test/energyLossFluctuation.png");

    // TOF.plotDeltaTime(Li6, "test/plotDeltaTime.png");

    // TOF.distributionOfReconstruction(Li6, 10000, true, true, "test/distributionOfReconstruction_linearMethod.png");
    // TOF.distributionOfReconstruction(Li6, 10000, false, true, "test/distributionOfReconstruction_nonLinearMethod.png");

    // TOF.plotDeltaBetaReciprocal(Li6, 0.4, 0.9, 5, true, "test/plotDeltaBetaReciprocal_linearMethod.png");
    TOF.plotDeltaBetaReciprocal(Li6, 0.4, 0.9, 5, false, "test/plotDeltaBetaReciprocal_nonLinearMethod.png");

#if configEnableDebug
    printf("[Info] The real 1 / beta: %f\n", 1 / Li6.getBeta());
#endif

    return 0;
}
