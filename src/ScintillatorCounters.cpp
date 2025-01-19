#include "ScintillatorCounters.h"

#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TLegend.h>

ScintillatorCounters::ScintillatorCounters(const double location, const bool direction, const double thickness, const double timeResolution, const Material &material)
    : location(location), direction(direction), thickness(thickness), timeResolution(timeResolution), material(material) {}

ScintillatorCounters::~ScintillatorCounters() {}

double ScintillatorCounters::energyLoss(const Particle &particle) const
{
    return this->material.linearStoppingPower(particle) * this->thickness;
}

void ScintillatorCounters::plotEnergyLoss(Particle particle, const double betaGammaMin, const double betaGammaMax, const int nPoints, const bool enableKineticEnergy, const std::string &fileName) const
{
    TCanvas *ScintillatorCounters_plotEnergyLossCanvas = new TCanvas("ScintillatorCounters_plotEnergyLossCanvas", "", 3508, 2480); // A4 size in pixels(300 dpi)
    ScintillatorCounters_plotEnergyLossCanvas->SetGrid();
    ScintillatorCounters_plotEnergyLossCanvas->cd();

    TMultiGraph *mg = new TMultiGraph();

    TGraph *graphEnergyLoss = new TGraph(nPoints + 1);
    graphEnergyLoss->SetLineColor(kBlue);

    TGraph *graphKineticEnergy = new TGraph(nPoints + 1);
    graphKineticEnergy->SetLineColor(kRed);

    const double minExponent = TMath::Log10(betaGammaMin);
    const double maxExponent = TMath::Log10(betaGammaMax);
    const double step = (maxExponent - minExponent) / nPoints;

    double yMax = 0;
    for (int i = 0; i <= nPoints; ++i)
    {
        const double betaGamma = TMath::Power(10, minExponent + i * step);
        particle.setBetaGamma(betaGamma);
        const double energyLoss = this->energyLoss(particle);
        yMax = TMath::Max(yMax, energyLoss);

        const double beta = particle.getBeta();
        graphEnergyLoss->SetPoint(i, beta, energyLoss);
        if (enableKineticEnergy)
            graphKineticEnergy->SetPoint(i, beta, particle.getEnergy() - particle.getMass0());
    }

    mg->Add(graphEnergyLoss);
    if (enableKineticEnergy)
        mg->Add(graphKineticEnergy);

    if (enableKineticEnergy)
        mg->SetTitle(";#beta;#DeltaE or E_{k} [MeV]");
    else
        mg->SetTitle(";#beta;#DeltaE [MeV]");

    mg->GetYaxis()->SetRangeUser(0, yMax);

    // QUESTION: When I put the definition of TCanvas here, the x-axis of the plot is not form betaGammaMin to betaGammaMax. Why?

    mg->Draw("AL");

    TLegend *legend = new TLegend(0.65, 0.75, 0.85, 0.85);
    legend->SetBorderSize(kNone);
    if (enableKineticEnergy)
    {
        legend->AddEntry(graphEnergyLoss, "#DeltaE (Energy loss)", "l");
        legend->AddEntry(graphKineticEnergy, "E_{k} (Kinetic energy)", "l");
        legend->Draw();
    }
    ScintillatorCounters_plotEnergyLossCanvas->SaveAs(fileName.c_str());

    return;
}

double ScintillatorCounters::LandauMostProbableEnergyLoss_xi(const Particle &particle) const
{
    constexpr double K = 0.307075;                                  // MeV mol^-1 cm^2 (4 * pi * N_A * r_e^2 * m_e * c^2, coefficient for dE/dx)
    const double z = particle.getCharge();                          // charge number of the particle (since we use charge in units of e, the charge number is the same as the charge)
    const double x = this->thickness * this->material.getDensity(); // x in g/cm^2
    const double Z = this->material.getZ();
    const double A = this->material.getA();
    const double beta = particle.getBeta();

    return (K / 2) * (Z / A) * (z * z) * (x / (beta * beta));
}

double ScintillatorCounters::LandauMostProbableEnergyLoss(const Particle &particle) const
{
    constexpr double massElectron = 0.51099895000; // Rest mass of the electron in MeV/c^2
    constexpr double j = 0.200;                    // Constant for the Landau most probable energy loss
    const double I = this->material.getI() * 1e-6; // convert eV to MeV
    const double beta = particle.getBeta();
    const double gamma = particle.getGamma();

    const double xi = this->LandauMostProbableEnergyLoss_xi(particle);

    const double partA = 2 * massElectron * beta * beta * gamma * gamma / I;
    const double partB = xi / I;

    return xi * (TMath::Log(partA) + TMath::Log(partB) + j - beta * beta - this->material.delta(beta, gamma));
}
