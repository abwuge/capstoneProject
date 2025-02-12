#include "ScintillatorCounters.h"

#include <TAxis.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>

ScintillatorCounters::ScintillatorCounters(
    const double    location,
    const bool      direction,
    const double    thickness,
    const double    timeResolution,
    const Material &material
)
    : sLocation(location), sDirection(direction), sThickness(thickness), sTimeResolution(timeResolution),
      sMaterial(material) {}

ScintillatorCounters::~ScintillatorCounters() {}

double ScintillatorCounters::energyLoss(const Particle &particle) const {
  return this->sMaterial.linearStoppingPower(particle) * this->sThickness;
}

TMultiGraph *ScintillatorCounters::energyLoss(
    Particle           particle,
    const double       betaGammaMin,
    const double       betaGammaMax,
    const int          nPoints,
    const bool         enableKineticEnergy,
    const bool         enablePlot,
    const std::string &fileName
) const {
  TCanvas *ScintillatorCounters_plotEnergyLossCanvas;
  if (enablePlot) {
    ScintillatorCounters_plotEnergyLossCanvas =
        new TCanvas("ScintillatorCounters_plotEnergyLossCanvas", "", 3508, 2480); // A4 size in pixels(300 dpi)
    ScintillatorCounters_plotEnergyLossCanvas->SetGrid();
    ScintillatorCounters_plotEnergyLossCanvas->cd();
  }

  TMultiGraph *mg = new TMultiGraph();

  TGraph *graphEnergyLoss = new TGraph(nPoints + 1);
  graphEnergyLoss->SetLineColor(kBlue);

  TGraph *graphKineticEnergy = nullptr;
  if (enableKineticEnergy) {
    graphKineticEnergy = new TGraph(nPoints + 1);
    graphKineticEnergy->SetLineColor(kRed);
  }

  const double minExponent = TMath::Log10(betaGammaMin);
  const double maxExponent = TMath::Log10(betaGammaMax);
  const double step        = (maxExponent - minExponent) / nPoints;

  double yMax = 0;
  for (int i = 0; i <= nPoints; ++i) {
    const double betaGamma = TMath::Power(10, minExponent + i * step);
    particle.setBetaGamma(betaGamma);

    double energyLoss;
    if (Config::useBetheBloch) energyLoss = this->energyLoss(particle);
    else
      energyLoss = this->LandauMostProbableEnergyLoss(particle);

    yMax = TMath::Max(yMax, energyLoss);

    const double beta = particle.getBeta();
    graphEnergyLoss->SetPoint(i, beta, energyLoss);
    if (enableKineticEnergy) graphKineticEnergy->SetPoint(i, beta, particle.getEnergy() - particle.getMass0());
  }

  mg->Add(graphEnergyLoss);
  if (enableKineticEnergy) mg->Add(graphKineticEnergy);

  if (enablePlot) {
    if (enableKineticEnergy) mg->SetTitle(";#beta;#DeltaE or E_{k} [MeV]");
    else
      mg->SetTitle(";#beta;#DeltaE [MeV]");

    mg->GetYaxis()->SetRangeUser(0, yMax);

    mg->Draw("AL");

    TLegend *legend = new TLegend(0.65, 0.75, 0.85, 0.85);
    legend->SetBorderSize(kNone);
    if (enableKineticEnergy) {
      legend->AddEntry(graphEnergyLoss, "#DeltaE (Energy loss)", "l");
      legend->AddEntry(graphKineticEnergy, "E_{k} (Kinetic energy)", "l");
      legend->Draw();
    }
    ScintillatorCounters_plotEnergyLossCanvas->SaveAs(fileName.c_str());

    delete legend;
    delete ScintillatorCounters_plotEnergyLossCanvas;
  }

  return mg;
}

double ScintillatorCounters::LandauMostProbableEnergyLoss_xi(const Particle &particle) const {
  double beta = particle.getBeta();
  if (beta < 1e-10) return 0;
  constexpr double K = 0.307075; // MeV mol^-1 cm^2 (4 * pi * N_A * r_e^2 * m_e * c^2, coefficient for dE/dx)
  const int        z = particle.getCharge(
  ); // charge number of the particle (since we use charge in units of e, the charge number is the same as the charge)
  const double x = this->sThickness * this->sMaterial.getDensity(); // x in g/cm^2
  const double Z = this->sMaterial.getZ();
  const double A = this->sMaterial.getA();

  return (K / 2) * (Z / A) * (z * z) * (x / (beta * beta));
}

double ScintillatorCounters::LandauMostProbableEnergyLoss(const Particle &particle) const {
  const double xi = this->LandauMostProbableEnergyLoss_xi(particle);
  if (!xi) return 0;

  constexpr double massElectron = 0.51099895000;                 // Rest mass of the electron in MeV/c^2
  constexpr double j            = 0.200;                         // Constant for the Landau most probable energy loss
  const double     I            = this->sMaterial.getI() * 1e-6; // convert eV to MeV
  const double     beta         = particle.getBeta();
  const double     gamma        = particle.getGamma();

  const double partA = 2 * massElectron * beta * beta * gamma * gamma / I;
  const double partB = xi / I;

  return xi * (TMath::Log(partA) + TMath::Log(partB) + j - beta * beta - this->sMaterial.delta(beta, gamma));
}

double ScintillatorCounters::LandauMostProbableEnergyLoss(const double xi, const Particle &particle) const {
  if (!xi) return 0;

  constexpr double massElectron = 0.51099895000;                 // Rest mass of the electron in MeV/c^2
  constexpr double j            = 0.200;                         // Constant for the Landau most probable energy loss
  const double     I            = this->sMaterial.getI() * 1e-6; // convert eV to MeV
  const double     beta         = particle.getBeta();
  const double     gamma        = particle.getGamma();

  const double partA = 2 * massElectron * beta * beta * gamma * gamma / I;
  const double partB = xi / I;

  return xi * (TMath::Log(partA) + TMath::Log(partB) + j - beta * beta - this->sMaterial.delta(beta, gamma));
}

void ScintillatorCounters::plotEnergyLossFluctuation(
    Particle           particle,
    const bool         enableKineticEnergy,
    TPad              *pad,
    const std::string &fileName
) const {
  const double xi                     = this->LandauMostProbableEnergyLoss_xi(particle);
  const double mostProbableEnergyLoss = this->LandauMostProbableEnergyLoss(xi, particle);
  const double sigma                  = 4.018 * xi;

  const double kineticEnergy = particle.getEnergy() - particle.getMass0();

  TF1 *landau = nullptr;
  if (xi) {
    if (enableKineticEnergy)
      landau = new TF1(
          "landau",
          "TMath::Landau(x, [0], [1])",
          mostProbableEnergyLoss - 3.5 * sigma,
          TMath::Max(mostProbableEnergyLoss + 20 * sigma, kineticEnergy + sigma)
      );
    else
      landau = new TF1(
          "landau",
          "TMath::Landau(x, [0], [1])",
          mostProbableEnergyLoss - 3.5 * sigma,
          mostProbableEnergyLoss + 20 * sigma
      );
    landau->SetParameters(mostProbableEnergyLoss, sigma);
  } else if (enableKineticEnergy) {
    landau = new TF1("landau", "[0]", -100, kineticEnergy + 100);
    landau->SetParameter(0, 0);
  } else {
    if (Config::enableWarning) {
      printf("[Warning] The most probable energy loss is 0! So, the Landau distribution is not plotted!\n");
      return;
    }
  }

  TCanvas *ScintillatorCounters_plotEnergyLossFluctuationCanvas;
  if (!pad) {
    ScintillatorCounters_plotEnergyLossFluctuationCanvas = new TCanvas(
        "ScintillatorCounters_plotEnergyLossFluctuationCanvas",
        "",
        3508,
        2480
    ); // A4 size in pixels(300 dpi)
    ScintillatorCounters_plotEnergyLossFluctuationCanvas->SetGrid();
    ScintillatorCounters_plotEnergyLossFluctuationCanvas->cd();
  } else {
    pad->SetGrid();
    pad->cd();
  }

  landau->SetTitle(";#DeltaE [MeV];");
  landau->Draw();

  TLegend *legend = new TLegend(0.65, 0.7, 0.85, 0.85);
  legend->SetBorderSize(kNone);
  legend->AddEntry("", Form("#Delta_{p} = %.4g", mostProbableEnergyLoss), "");
  legend->AddEntry("", Form("#sigma = %.4g", sigma), "");
  TLine *kineticEnergyLine =
      new TLine(particle.getEnergy() - particle.getMass0(), 0, particle.getEnergy() - particle.getMass0(), 1);
  if (enableKineticEnergy) {
    legend->AddEntry(kineticEnergyLine, "Kinetic energy", "l");
    kineticEnergyLine->Draw("same");
  }

  legend->Draw();

  if (!pad) {
    ScintillatorCounters_plotEnergyLossFluctuationCanvas->SaveAs(fileName.c_str());
    delete kineticEnergyLine;
    delete legend;
    delete landau;
    delete ScintillatorCounters_plotEnergyLossFluctuationCanvas;
  }
}
