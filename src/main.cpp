#include <algorithm>
#include <iostream>
#include <thread>
#include <tuple>
#include <vector>

#include <TAxis.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TROOT.h>
#include <TVector3.h>

#include "Config.h"
#include "Detector.h"
#include "Material.h"
#include "Particle.h"
#include "ScintillatorCounters.h"
#include "ThreadPool.h"

int main(int argc, char *argv[]) {
  /* BEGIN Initialization */
  gROOT->SetStyle("Pub");

  Config::enableFixedSeed = true;
  Config::config          = Config::getInstance();
  Config::random          = Config::getRandom();

  if (Config::enableMultiThreading) {
    ROOT::EnableThreadSafety();

    if (argc > 1)
      ThreadPool::getInstance(std::min(std::thread::hardware_concurrency(), (unsigned int)std::stoi(argv[1])));
    else
      ThreadPool::getInstance(std::thread::hardware_concurrency());

    printf("[Info] Number of threads: %lu\n", ThreadPool::getInstance().getNumThreads());
  }
  /* END Initialization */

  /* BEGIN Material properties */
  const Material *polyvinyltoluene;
  if (Config::todo) // IDK some porperties of polyvinyltoluene, so I used polystyrene to test the code
    polyvinyltoluene = new Material(materials.at(MaterialName::Polystyrene));
  else
    // The EJ-200 scintillator counter is made of polyvinyltoluene [(2-CH3C6H4CHCH2)n]
    polyvinyltoluene = new Material(materials.at(MaterialName::Polyvinyltoluene)); // Polyvinyltoluene material
  /* END Material properties */

  /* BEGIN Detector properties */
  // z-coordinate of EJ-200 scintillator counters in AMS-02 in cm
  // (assuming they are infinite planes in the x and y directions)
  const std::vector<double> EJ_200Location = {-65, -63, 63, 65};

  // Direction of EJ-200 scintillator counters in AMS-02 (1 for x, 0 for y)
  const std::vector<bool> EJ_200Direction = {1, 0, 0, 1};
  constexpr double        thickness       = 1;   // Thickness of EJ-200 scintillator counters in cm
  constexpr double        timeResolution  = 0.1; // Time resolution of EJ-200 scintillator counters in ns

  std::vector<ScintillatorCounters> EJ_200;      // EJ-200 scintillator counters in AMS-02
  EJ_200.reserve(4);
  for (int i = 0; i < 4; i++) {
    EJ_200.push_back(
        ScintillatorCounters(EJ_200Location[i], EJ_200Direction[i], thickness, timeResolution, *polyvinyltoluene)
    );
  }
  // const TVector3 B(0, 0.14, 0); // Magnetic field in AMS-02 in Tesla
  const TVector3 B(0, 0, 0);     // used for testing

  const Detector TOF(EJ_200, B); // Time-of-flight detector in AMS-02
  /* END Detector properties */

  constexpr double uIn_kg = 1.66053906660e-27; // Atomic mass unit in kg
                                               // (from https://www.bipm.org/documents/20126/41483022/SI-Brochure-9.pdf)
  constexpr double kg = TMath::C() * TMath::C() / (1e6 * TMath::Qe()); // 1 kg in MeV/c^2
  constexpr double u  = uIn_kg * kg;                                   // Atomic mass unit in MeV/c^2
  /* BEGIN Particle Li6 properties */
  // The initial position of the particle is fixed at the center of the first scintillator counter
  constexpr double chargeLi6    = 3;              // Charge of Li6 in e
  constexpr double mass0Li6     = 6.01512289 * u; // Rest mass of Li6 in MeV/c^2 (from https://ciaaw.org/lithium.htm)
  constexpr double startBetaLi6 = 0.4;            // Initial beta of Li6
  const TVector3   startPositionLi6(0, 0, TOF.getMinZ());            // Initial position of Li6 in cm

  Particle Li6(chargeLi6, mass0Li6, startBetaLi6, startPositionLi6); // Initial Li6
  /* END Particle Li6 properties */

  /* BEGIN Particle Proton properties */
  // The initial position of the particle is fixed at the center of the first scintillator counter
  constexpr double chargeProton = 1;      // Charge of proton in e
  constexpr double mass0Proton =
      1.67262192595e-27 * kg;             // Rest mass of proton in MeV/c^2
                                          //(from https://physics.nist.gov/cgi-bin/cuu/Value?mp|search_for=proton+mass)
  constexpr double startBetaProton = 0.4; // Initial beta of proton
  const TVector3   startPositionProton(0, 0, TOF.getMinZ());                        // Initial position of proton in cm

  Particle proton(chargeProton, mass0Proton, startBetaProton, startPositionProton); // Initial proton
  /* END Particle Proton properties */

  // Draw the energy loss of the particle in the scintillator counters using Bethe-Bloch and Landau
  if (false) {
    TCanvas *canvas = new TCanvas("canvas", "", 3508, 2480); // A4 size in pixels(300 dpi)
    canvas->SetGrid();
    canvas->cd();

    TLegend *legend = new TLegend(0.65, 0.75, 0.85, 0.85);
    legend->SetBorderSize(kNone);

    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle(";#beta;#DeltaE or E_{k} [MeV]");

    Config::useBetheBloch             = true;
    TMultiGraph *multiGraphBetheBloch = EJ_200.at(0).energyLoss(Li6, 0.1, 1000, 1000, true);
    TGraph      *graphBetheBloch      = (TGraph *)multiGraphBetheBloch->GetListOfGraphs()->At(0);
    graphBetheBloch->SetLineColor(kBlue);
    mg->Add(graphBetheBloch);
    legend->AddEntry(graphBetheBloch, "Bethe-Bloch", "l");
    TGraph *graphKineticEnergy = (TGraph *)multiGraphBetheBloch->GetListOfGraphs()->At(1);
    graphKineticEnergy->SetLineColor(kRed);
    mg->Add(graphKineticEnergy);
    legend->AddEntry(graphKineticEnergy, "Kinetic Energy", "l");

    Config::useBetheBloch         = false;
    TMultiGraph *multiGraphLandau = EJ_200.at(0).energyLoss(Li6, 0.1, 1000, 1000, false);
    TGraph      *graphLandau      = (TGraph *)multiGraphLandau->GetListOfGraphs()->At(0);
    graphLandau->SetLineColor(kGreen);
    mg->Add(graphLandau);
    legend->AddEntry(graphLandau, "Landau", "l");

    double yMax = 0;
    for (int i = 0; i < graphBetheBloch->GetN(); ++i) yMax = TMath::Max(yMax, graphBetheBloch->GetY()[i]);
    for (int i = 0; i < graphLandau->GetN(); ++i) yMax = TMath::Max(yMax, graphLandau->GetY()[i]);

    mg->GetYaxis()->SetRangeUser(0, yMax);

    mg->Draw("AL");
    legend->Draw();
    canvas->SaveAs("test/energyLoss.png");

    mg->GetListOfGraphs()->Clear();
    delete mg;
    delete multiGraphBetheBloch;
    delete multiGraphLandau;
    delete canvas;
    delete legend;
  }

  // Draw the energy loss distribution of the particle in the scintillator counters using Landau
  if (false) {
    TCanvas *canvas = new TCanvas("canvas", "", 3508, 2480); // A4 size in pixels(300 dpi)
    canvas->Divide(2, 2);

    if (Config::enableDebug) printf("[Info] Li6 energy: %f MeV, velocity: %f c\n", Li6.getEnergy(), Li6.getBeta());
    for (int i = 0; i < 4; ++i) {
      EJ_200.at(0).plotEnergyLossFluctuation(Li6, 100000, true, (TPad *)canvas->GetPad(i + 1));
      const double xi = EJ_200.at(i).LandauMostProbableEnergyLoss_xi(Li6);
      const double el = EJ_200.at(i).LandauMostProbableEnergyLoss(xi, Li6);
      if (Config::enableEnergyLossFluctuation) Li6.setEnergy(Li6.getEnergy() - Config::random->Landau(el, 4.018 * xi));
      else
        Li6.setEnergy(Li6.getEnergy() - el);
      if (Config::enableDebug) printf("[Info] Li6 energy: %f MeV, velocity: %f c\n", Li6.getEnergy(), Li6.getBeta());
    }

    canvas->SaveAs("test/energyLossFluctuation.png");

    delete canvas;

    Li6.setBeta(startBetaLi6);
  }

  // Draw the kinetic energy distribution after random Landau energy loss
  if (false) {
    TCanvas *canvas = new TCanvas("canvas", "", 3508, 2480); // A4 size in pixels(300 dpi)
    canvas->Divide(2, 2);

    std::vector<TH1F *> histograms(4);
    for (int i = 0; i < 4; ++i)
      histograms.at(i) = new TH1F(Form("TOF%d", i), ";E_{k} [MeV];", 100, 0, Li6.getEnergy() - Li6.getMass0() + 1);

    std::vector<TLegend *> legends(4);
    std::vector<int>       cnt(4, 0);

    auto Add = [&](Particle particle) {
      for (int i = 0; i < 4; ++i) {
        const double xi = EJ_200.at(i).LandauMostProbableEnergyLoss_xi(particle);
        const double el = EJ_200.at(i).LandauMostProbableEnergyLoss(xi, particle);
        if (Config::enableEnergyLossFluctuation)
          particle.setEnergy(particle.getEnergy() - Config::random->Landau(el, 4.018 * xi));
        else
          particle.setEnergy(particle.getEnergy() - el);

        double energy = particle.getEnergy() - particle.getMass0();
        histograms.at(i)->Fill(energy);
        if (energy == 0) break;
        cnt[i]++;
      }
    };

    for (int i = 0; i < 10000; ++i) Add(Li6);

    for (int i = 0; i < 4; ++i) {
      canvas->cd(i + 1);
      histograms.at(i)->Draw();

      legends.at(i) = new TLegend(0.18, 0.78, 0.38, 0.87);
      legends.at(i)->SetBorderSize(kNone);
      legends.at(i)->AddEntry("", Form("Entries: %d", (int)histograms.at(i)->GetEntries()), "");
      legends.at(i)->AddEntry("", Form("E_{k} > 0: %d", cnt.at(i)), "");
      legends.at(i)->Draw();
    }

    canvas->SaveAs("test/kineticEnergyDistribution.png");

    for (auto legend : legends) delete legend;
    for (auto histogram : histograms) delete histogram;
    delete canvas;
  }

  // TOF.plotDeltaTime(Li6, "test/plotDeltaTime.png");

  // TOF.distributionOfReconstruction(Li6, 10000, true, true, "test/distributionOfReconstruction_linearMethod.png");
  if (true) {
    proton.setBeta(0.4);
    TOF.distributionOfReconstruction(
        proton,
        10000,
        false,
        true,
        "test/distributionOfReconstruction_nonLinearMethod.png"
    );
  }

  // Draw the difference between real and reconstructed 1/beta of the particle in this detectors
  if (false) {
    TCanvas *canvas = new TCanvas("canvas", "", 3508, 2480); // A4 size in pixels(300 dpi)
    canvas->SetGrid();
    canvas->cd();

    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle(";#beta;#Delta(1/#beta)");

    TLegend *legend = new TLegend(0.45, 0.2, 0.88, 0.32);
    legend->SetBorderSize(kNone);

    const int       nPoints        = 5;
    const Particle &detectParticle = proton;

    Config::useLandau = true;
    // Linear method (Landau Fluctuation)
    if (false) {
      TGraphErrors *graphErrorsLinearMethodLandauFluctuation =
          TOF.deltaBetaReciprocal(detectParticle, 0.4, 0.9, nPoints, true);
      graphErrorsLinearMethodLandauFluctuation->SetLineColor(kBlue);
      mg->Add(graphErrorsLinearMethodLandauFluctuation);
      legend->AddEntry(graphErrorsLinearMethodLandauFluctuation, "Linear Method (Landau Fluctuation)", "l");
    }

    // Non-linear method (Landau Fluctuation)
    if (true) {
      TGraphErrors *graphErrorsNonLinearMethodLandauFluctuation =
          TOF.deltaBetaReciprocal(detectParticle, 0.4, 0.9, nPoints, false);
      graphErrorsNonLinearMethodLandauFluctuation->SetLineColor(kRed);
      mg->Add(graphErrorsNonLinearMethodLandauFluctuation);
      legend->AddEntry(graphErrorsNonLinearMethodLandauFluctuation, "Non-Linear Method (Landau Fluctuation)", "l");
    }

    Config::useLandau = false;
    // Non-linear method (Gaussian Fluctuation)
    if (true) {
      TGraphErrors *graphErrorsNonLinearMethodGaussianFluctuation =
          TOF.deltaBetaReciprocal(detectParticle, 0.4, 0.9, nPoints, false);
      graphErrorsNonLinearMethodGaussianFluctuation->SetLineColor(kGreen);
      mg->Add(graphErrorsNonLinearMethodGaussianFluctuation);
      legend->AddEntry(graphErrorsNonLinearMethodGaussianFluctuation, "Non-Linear Method (Gaussian Fluctuation)", "l");
    }

    mg->Draw("AL");
    legend->Draw();

    canvas->SaveAs("test/deltaBetaReciprocal.png");

    mg->GetYaxis()->SetRangeUser(-0.12, 0.1);
    legend->SetY1NDC(0.76);
    legend->SetY2NDC(0.88);

    canvas->SaveAs("test/deltaBetaReciprocal_zoom.png");

    delete canvas;
    delete mg;
    delete legend;
  }

  if (Config::enableDebug) printf("[Info] The real 1 / beta: %f\n", 1 / Li6.getBeta());

  return 0;
}
