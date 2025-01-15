#include <iostream>
#include <tuple>
#include <vector>

#include <TApplication.h>
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

#if configTodo // IDK some porperties of polyvinyltoluene, so I used polystyrene to test the code
    /* BEGIN Material properties */
    const Material polyvinyltoluene(materials.at(MaterialName::Polystyrene));
    /* END Material properties */
#else
    /* BEGIN Material properties */
    // The EJ-200 scintillator counter is made of polyvinyltoluene [(2-CH3C6H4CHCH2)n]
    const Material polyvinyltoluene(materials.at(MaterialName::Polyvinyltoluene)); // Polyvinyltoluene material
    /* END Material properties */
#endif

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
    // const TVector3 B(0, 0.14, 0); // Magnetic field in AMS-02 in Tesla
    const TVector3 B(0, 0, 0); // used for testing

    const Detector TOF(EJ_200, B); // Time-of-flight detector in AMS-02
    /* END Detector properties */

    /* BEGIN Particle properties */
    // The initial position of the particle is fixed at the center of the first scintillator counter
    constexpr double charge = 3;                       // Charge of Li6 in e
    constexpr double u = 931.49410372;                 // Atomic mass unit in MeV/c^2 (from https://en.wikipedia.org/wiki/Dalton_(unit))
    constexpr double mass0 = 6.0151228874 * u;         // Rest mass of Li6 in MeV/c^2 (from https://en.wikipedia.org/wiki/Isotopes_of_lithium)
    constexpr double startBeta = 0.3;                  // Initial beta of Li6
    const TVector3 startPosition(0, 0, TOF.getMinZ()); // Initial position of Li6 in cm

    Particle Li6(charge, mass0, startBeta, startPosition); // Initial Li6
    Li6.setBetaGamma(10);                                  // Set the beta * gamma of Li6
    std::cout << EJ_200.at(0).energyLoss(Li6) << std::endl;
    /* END Particle properties */

    std::vector<std::tuple<double, double, TVector3>> hitDataWithEnergyLoss = TOF.particleHitData(Li6, true);
    std::vector<std::tuple<double, double, TVector3>> hitDataWithoutEnergyLoss = TOF.particleHitData(Li6, false);

    // std::cout << "Hit data with energy loss:" << std::endl;
    // for (const auto &hit : hitDataWithEnergyLoss)
    //     std::cout << "Hit time: " << std::get<0>(hit) << " ns, Propagation length: " << std::get<1>(hit) << " cm, Hit position: (" << std::get<2>(hit).X() << ", " << std::get<2>(hit).Y() << ", " << std::get<2>(hit).Z() << ") cm" << std::endl;
    // std::cout << "Hit data without energy loss:" << std::endl;
    // for (const auto &hit : hitDataWithoutEnergyLoss)
    //     std::cout << "Hit time: " << std::get<0>(hit) << " ns, Propagation length: " << std::get<1>(hit) << " cm, Hit position: (" << std::get<2>(hit).X() << ", " << std::get<2>(hit).Y() << ", " << std::get<2>(hit).Z() << ") cm" << std::endl;

    EJ_200.at(0).plotEnergyLoss(Li6);

    TGraph *graph = new TGraph();
    graph->SetTitle("#Deltat vs. Propagation length;Propagation length [cm];#Deltat [ns]");

    for (int i = 0; i < TOF.getScintillatorCounters().size(); ++i)
        try
        {
            graph->SetPoint(i, std::get<1>(hitDataWithEnergyLoss.at(i)), std::get<0>(hitDataWithEnergyLoss.at(i)) - std::get<0>(hitDataWithoutEnergyLoss.at(i)));
        }
        catch (const std::out_of_range &e)
        {
#if configEnableWarning
            printf("[Warning] Seems like the particle has not hit the scintillator counter since counter %d! Only plotting the hit scintillator counters!\n", i);
#endif
            break;
        }

    TCanvas *c1 = new TCanvas("c1InMain", "", 3508, 2480); // A4 size in pixels(300 dpi)
    c1->cd();

    graph->Draw("AL");

    c1->SaveAs("test2.png");

    // std::cout << "Cyclotron radius: " << TOF.particleCyclotronRadius(Li6) << " cm" << std::endl;
    // TVector3 direction = TOF.particleCyclotronDirection(Li6);
    // std::cout << "Cyclotron direction: (" << direction.X() << ", " << direction.Y() << ", " << direction.Z() << ")" << std::endl;
    // std::cout << "beta: " << Li6.getBeta() << std::endl;
    // std::cout << "Stopping power: " << EJ_200.at(0).energyLoss(Li6) << " MeV" << std::endl;
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
