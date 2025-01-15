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
    TCanvas *c1 = new TCanvas("c1InScintillatorCounters::plotEnergyLoss", "", 3508, 2480); // A4 size in pixels(300 dpi)
    c1->SetLogx();

    TMultiGraph *mg = new TMultiGraph();

    TGraph *graphEnergyLoss = new TGraph(nPoints + 1);
    graphEnergyLoss->SetLineColor(kBlue);
    graphEnergyLoss->SetLineWidth(3);

    TGraph *graphKineticEnergy = new TGraph(nPoints);
    graphKineticEnergy->SetLineColor(kRed);
    graphKineticEnergy->SetLineWidth(3);

    const double minExponent = TMath::Log10(betaGammaMin);
    const double maxExponent = TMath::Log10(betaGammaMax);
    const double step = (maxExponent - minExponent) / nPoints;

    double yMax = 0;
    for (int i = 0; i <= nPoints; ++i)
    {
        double betaGamma = TMath::Power(10, minExponent + i * step);
        particle.setBetaGamma(betaGamma);
        double energyLoss = this->energyLoss(particle);
        yMax = TMath::Max(yMax, energyLoss);
        graphEnergyLoss->SetPoint(i, betaGamma, energyLoss);
        if (enableKineticEnergy)
            graphKineticEnergy->SetPoint(i, betaGamma, particle.getEnergy() - particle.getMass0());
    }

    mg->Add(graphEnergyLoss);
    if (enableKineticEnergy)
        mg->Add(graphKineticEnergy);

    if (enableKineticEnergy)
        mg->SetTitle(";#beta#gamma;Energy loss / Kinetic energy [MeV]");
    else
        mg->SetTitle(";#beta#gamma;Energy loss [MeV]");

    mg->GetYaxis()->SetRangeUser(0, yMax);

    // QUESTION: When I put the definition of c1 here, the x-axis of the plot is not form betaGammaMin to betaGammaMax. Why?
    // TCanvas *c1 = new TCanvas("c1InScintillatorCounters::plotEnergyLoss", "", 3508, 2480); // A4 size in pixels(300 dpi)
    // c1->SetLogx();

    c1->cd();

    mg->Draw("AL");

    TLegend *legend = new TLegend(0.7, 0.8, 0.9, 0.9);
    if (enableKineticEnergy)
    {
        legend->AddEntry(graphEnergyLoss, "Energy loss", "l");
        legend->AddEntry(graphKineticEnergy, "Kinetic energy", "l");
        legend->Draw();
    }
    c1->SaveAs(fileName.c_str());

    return;
}
