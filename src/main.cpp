#include <iostream>
#include <tuple>
#include <vector>

#include <TROOT.h>
#include <TVector3.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TAxis.h>
#include <TLegend.h>

#include "config.h"
#include "Material.h"
#include "ScintillatorCounters.h"
#include "Particle.h"
#include "Detector.h"

void test(int argc, char *argv[]);

int main(int argc, char *argv[])
{
#if configEnableTest
    test(argc, argv);
    return 0;
#endif

    gROOT->SetStyle("Pub");

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
    /* END Detector properties */

    /* BEGIN Particle properties */
    // The initial position of the particle is fixed at the center of the first scintillator counter
    constexpr double charge = 3;                                                 // Charge of Li6 in e
    constexpr double uIn_kg = 1.66053906660e-27;                                 // Atomic mass unit in kg (from https://www.bipm.org/documents/20126/41483022/SI-Brochure-9.pdf)
    constexpr double u = uIn_kg * TMath::C() * TMath::C() / (1e6 * TMath::Qe()); // Atomic mass unit in MeV/c^2
    constexpr double mass0 = 6.01512289 * u;                                     // Rest mass of Li6 in MeV/c^2 (from https://ciaaw.org/lithium.htm)
    constexpr double startBeta = 0.7;                                            // Initial beta of Li6
    const TVector3 startPosition(0, 0, TOF.getMinZ());                           // Initial position of Li6 in cm

    Particle Li6(charge, mass0, startBeta, startPosition); // Initial Li6
    /* END Particle properties */

    // plot the energy loss of Li6 in EJ-200 scintillator counters
    EJ_200.at(0).plotEnergyLoss(Li6, 0.1, 1000, 1000, true, "result/energyLoss.png");

    std::vector<std::tuple<double, double, TVector3>> hitDataWithEnergyLoss = TOF.particleHitData(Li6, true);
    std::vector<std::tuple<double, double, TVector3>> hitDataWithoutEnergyLoss = TOF.particleHitData(Li6, false);

    TOF.plotDeltaTime(hitDataWithEnergyLoss, hitDataWithoutEnergyLoss, "result/deltaTime.png");

    // std::vector<double> hitTimes, propagationLengths;

    // hitTimes.reserve(hitDataWithEnergyLoss.size()), propagationLengths.reserve(hitDataWithEnergyLoss.size());
    // for (const auto &hitData : hitDataWithEnergyLoss)
    // {

    //     hitTimes.push_back(std::get<0>(hitData));
    //     propagationLengths.push_back(std::get<1>(hitData));
    // }

    TOF.plotReconstructDataUsingLinearMethod(Li6, "result/reconstructData.png");

    TOF.distributionOfReconstructionUsingLinearMethod(Li6, 10000, true, "result/distributionOfReconstruction.png");

    TOF.plotDeltaBetaReciprocal(Li6, 0.4, 0.9, 5, "result/deltaBetaReciprocal.png");

#if configEnableDebug
    printf("[Info] The real 1 / beta: %f\n", 1 / Li6.getBeta());
#endif

    return 0;
}

#if configEnableTest

void test(int argc, char *argv[])
{
    /* BEGIN Material properties */
    const Material copper(materials.at(MaterialName::Copper)); // Copper material
    /* END Material properties */

    /* BEGIN Material properties 2 */
    const Material polystyrene(materials.at(MaterialName::Polystyrene)); // Polystyrene material
    /* END Material properties 2 */

    /* BEGIN Particle properties */
    constexpr double chargeElectron = -1;            // Charge in e
    constexpr double mass0Electron = 0.5109989461;   // Rest mass in MeV/c^2
    const TVector3 startMomentumElectron(0, 0, 1e3); // Initial momentum in MeV/c
    const TVector3 startPositionElectron(0, 0, 0);   // Initial position in cm

    Particle electron(chargeElectron, mass0Electron, startMomentumElectron, startPositionElectron); // Positive muon
    /* END Particle properties */

    /* BEGIN Particle properties 2 */
    constexpr double chargeLi6 = 3;               // Charge in e
    constexpr double u = 931.49410372;            // Atomic mass unit in MeV/c^2 (from https://en.wikipedia.org/wiki/Dalton_(unit))
    constexpr double mass0Li6 = 6.0151228874 * u; // Rest mass in MeV/c^2
    const TVector3 startMomentumLi6(0, 0, 1e3);   // Initial momentum in MeV/c
    const TVector3 startPositionLi6(0, 0, 0);     // Initial position in cm

    Particle Li6(chargeLi6, mass0Li6, startMomentumLi6, startPositionLi6); // Positive muon
    /* END Particle properties 2 */

    TCanvas *c1 = new TCanvas("c1", "", 3508, 2480); // A4 size in pixels(300 dpi)
    TGraph *graph = new TGraph();
    graph->SetTitle("Li6 on Polystyrene;#beta#gamma;Mass Stopping Power [MeV cm^{2}/g]");
    c1->SetLogx();

    for (int i = 0; i <= 1000; ++i)
    {
        double betaGamma = TMath::Power(10, -1 + 4 * i / 1000.);
        Li6.setBetaGamma(betaGamma);
        graph->SetPoint(i, Li6.getBeta() * Li6.getGamma(), polystyrene.massStoppingPower(Li6));
    }

    graph->Draw("AL");
    c1->SaveAs("test.png");
}

#endif
