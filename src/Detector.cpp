#include "Detector.h"

#include <algorithm>

#include <TMath.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TF1.h>
#include <TH1F.h>

// Infact, we can use TRandom3 *Detector::Random = new TRandom3(configEnableFixedSeed); to implement the following code
// But the following code is more readable
#if configEnableFixedSeed
TRandom3 *Detector::Random = new TRandom3(1); // Fixed seed 1
#else
TRandom3 *Detector::Random = new TRandom3(0); // 0 means time-based seed
#endif

Detector::Detector(const std::vector<ScintillatorCounters> &scintillatorCounters, const TVector3 &B)
    : scintillatorCounters(scintillatorCounters), B(B)
{
    std::sort(this->scintillatorCounters.begin(), this->scintillatorCounters.end(),
              [](const ScintillatorCounters &a, const ScintillatorCounters &b)
              { return a.getLocation() < b.getLocation(); });
}

Detector::~Detector() {}

double Detector::particleCyclotronRadius(const Particle &particle) const
{
    double absCharge = std::abs(particle.getCharge());
    const TVector3 &momentum = particle.getMomentum();

    double radius = momentum.Perp(this->B) / (absCharge * this->B.Mag()); // cyclotron radius in MeV / c / (e * T) = 1e6 J / (c * T)
    constexpr double conversionFactor = 1e6 / TMath::C() * 1e2;           // conversion factor from MeV / c / (e * T) to cm
    return radius * conversionFactor;
}

TVector3 Detector::particleCyclotronDirection(const Particle &particle) const
{
    double charge = particle.getCharge();
    double chargeSign = charge / std::abs(charge);
    const TVector3 &momentumUnit = particle.getMomentum().Unit();
    const TVector3 &BUint = this->B.Unit();

    return momentumUnit.Cross(BUint) * chargeSign;
}

std::vector<std::tuple<double, double, TVector3>> Detector::particleHitData(Particle particle, const bool enableEnergyLoss) const
{
    std::vector<std::tuple<double, double, TVector3>> hitData;
    hitData.reserve(this->scintillatorCounters.size());

    if (this->B.Mag() == 0)
    {
        constexpr double conversionFactor = TMath::Ccgs() * 1e-9; // conversion factor from c to cm/ns
        // The trajectory is a straight line (x, y, z) = (x0, y0, z0) + t * (vx, vy, vz)
        // The hit time is calculated by the formula: t = (z - z0) / vz
        double time = 0;
        double propagationLength = 0;
        for (const ScintillatorCounters &scintillatorCounter : this->scintillatorCounters)
        {
            const double z0 = particle.getPosition().Z();
            const double vz = particle.getVelocity().Z() * conversionFactor; // velocity in cm/ns
            if (vz <= 0)
            {
                printf("Cannot hit any more scintillator counters! The particle is moving in the opposite direction!\n");
                break;
            }
            const double z = scintillatorCounter.getLocation();
            const double deltaTime = (z - z0) / vz;
            const TVector3 deltaX = deltaTime * particle.getVelocity() * conversionFactor;
            const TVector3 hitPosition = particle.getPosition() + deltaX;

            particle.setPosition(hitPosition);
            if (enableEnergyLoss)
            {
                const double particleEnergyLoss = scintillatorCounter.energyLoss(particle);
                particle.setEnergy(particle.getEnergy() - particleEnergyLoss);
#if configEnableDebug
                printf("[Info] Energy loss: %f MeV, Energy: %f MeV, Velocity: %f c\n", particleEnergyLoss, particle.getEnergy(), particle.getVelocity().Mag());
#endif
            }

            time += deltaTime, propagationLength += deltaX.Mag();
            hitData.push_back({time, propagationLength, hitPosition});
        }
    }
    else
    {
        printf("[Error] Have not implemented the case when the magnetic field is not zero!\n");
        exit(1);
        // const double radius = particleCyclotronRadius(particle);
        // const TVector3 radiusVector = radius * particleCyclotronDirection(particle);
        // const TVector3 center = particle.getPosition() + radiusVector;

        // printf("radius: %f\n", radius);

        // const TVector3 &momentum = particle.getMomentum();
        // const TVector3 &position = particle.getPosition();
        // const TVector3 &velocity = particle.getVelocity();

        // for (const ScintillatorCounters &scintillatorCounter : scintillatorCounters)
        // {
        //     const double thickness = scintillatorCounter.getThickness();
        //     const Material &material = scintillatorCounter.getMaterial();

        //     const double beta = particle.getBeta();
        //     const double gamma = particle.getGamma();
        //     const double delta = material.delta(beta, gamma);
        //     const double massStoppingPower = material.massStoppingPower(particle);
        //     const double energyLoss = massStoppingPower * thickness * (1 - delta);

        //     const double hitTime = thickness / (velocity.Mag() * TMath::Cos(velocity.Angle(momentum)));
        //     const TVector3 hitPosition = position + hitTime * velocity;

        //     hitData.push_back(std::make_pair(hitTime, hitPosition));
        // }
    }

    return hitData;
}

void Detector::plotDeltaTime(const Particle &particle, const std::string &fileName) const
{
    std::vector<std::tuple<double, double, TVector3>> hitDataWithEnergyLoss = this->particleHitData(particle, true);
    std::vector<std::tuple<double, double, TVector3>> hitDataWithoutEnergyLoss = this->particleHitData(particle, false);

    this->plotDeltaTime(hitDataWithEnergyLoss, hitDataWithoutEnergyLoss, fileName);
}

void Detector::plotDeltaTime(const std::vector<std::tuple<double, double, TVector3>> &hitDataWithEnergyLoss, const std::vector<std::tuple<double, double, TVector3>> &hitDataWithoutEnergyLoss, const std::string &fileName) const
{
    TCanvas *Detector_plotDeltaTimeCanvas = new TCanvas("Detector_plotDeltaTimeCanvas", "", 3508, 2480); // A4 size in pixels(300 dpi)
    Detector_plotDeltaTimeCanvas->cd();

    TGraph *graph = new TGraph(this->scintillatorCounters.size());
    graph->SetTitle("#Deltat vs. Propagation length;Propagation length [cm];#Deltat [ns]");

    for (int i = 0; i < this->scintillatorCounters.size(); ++i)
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

    graph->Draw("AL");

    Detector_plotDeltaTimeCanvas->SaveAs(fileName.c_str());
}

std::vector<double> Detector::detect(const Particle &particle) const
{
    std::vector<std::tuple<double, double, TVector3>> hitData = this->particleHitData(particle, true);
    std::vector<double> detectedTimes;
    detectedTimes.reserve(hitData.size());

    for (int i = 0; i < hitData.size(); ++i)
    {
        const double hitTime = std::get<0>(hitData.at(i));
        detectedTimes.push_back(this->Random->Gaus(hitTime, this->scintillatorCounters.at(i).getTimeResolution()));
    }

    return detectedTimes;
}

std::vector<double> Detector::detect(const std::vector<double> &hitTimes) const
{
    std::vector<double> detectedTimes;
    detectedTimes.reserve(hitTimes.size());

    for (int i = 0; i < hitTimes.size(); ++i)
    {
        const double hitTime = hitTimes.at(i);
        detectedTimes.push_back(this->Random->Gaus(hitTime, this->scintillatorCounters.at(i).getTimeResolution()));
    }

    return detectedTimes;
}

double Detector::reconstructUsingLinearMethod(const Particle &particle) const
{
    const std::vector<std::tuple<double, double, TVector3>> hitData = this->particleHitData(particle, true);
    std::vector<double> hitTimes, propagationLengths;
    hitTimes.reserve(hitData.size()), propagationLengths.reserve(hitData.size());

    for (const auto &hit : hitData)
    {
        hitTimes.push_back(std::get<0>(hit));
        propagationLengths.push_back(std::get<1>(hit));
    }

    return this->reconstructUsingLinearMethod(this->detect(hitTimes), propagationLengths);
}

double Detector::reconstructUsingLinearMethod(const std::vector<double> &detectedTimes, const std::vector<double> &propagationLengths) const
{
    const int n = detectedTimes.size();
    if (n != propagationLengths.size())
    {
        printf("[Error] The sizes of detected times and propagation lengths are not the same!\n");
        exit(1);
    }

    TGraph *graph = new TGraph(n, &propagationLengths[0], &detectedTimes[0]);

    TF1 *f1 = new TF1("f1", "[0] * x", 0, propagationLengths.back());

    graph->Fit(f1, "Q");

    constexpr double conversionFactor = TMath::Ccgs() * 1e-9; // conversion factor from ns/cm to 1/c
    const double betaReciprocal = f1->GetParameter(0) * conversionFactor;

#if configEnableDebug
    printf("[Info] The reconstructed 1/beta using the linear method: %f\n", betaReciprocal);
#endif

    return betaReciprocal;
}

void Detector::plotReconstructDataUsingLinearMethod(const Particle &particle, const std::string &fileName) const
{
    const std::vector<std::tuple<double, double, TVector3>> hitData = this->particleHitData(particle, true);
    const int n = hitData.size();
    std::vector<double> hitTimes, propagationLengths;
    hitTimes.reserve(n), propagationLengths.reserve(n);

    for (const auto &hit : hitData)
    {
        hitTimes.push_back(std::get<0>(hit));
        propagationLengths.push_back(std::get<1>(hit));
    }

    std::vector<double> detectedTimes = this->detect(hitTimes);

    TGraph *graphRealData = new TGraph(n, &propagationLengths[0], &hitTimes[0]);
    graphRealData->SetLineWidth(3);

    TGraph *graphReconstructData = new TGraph(n, &propagationLengths[0], &detectedTimes[0]);
    graphReconstructData->SetMarkerStyle(20);
    graphReconstructData->SetMarkerSize(3);
    graphReconstructData->SetTitle(";Propagation length [cm];Time [ns]");

    TCanvas *Detector_plotReconstructDataCanvas = new TCanvas("Detector_plotReconstructDataCanvas", "", 3508, 2480); // A4 size in pixels(300 dpi)
    Detector_plotReconstructDataCanvas->cd();

    TF1 *f1 = new TF1("f1", "[0] * x", 0, propagationLengths.back());

    graphReconstructData->Draw("AP");
    graphReconstructData->Fit(f1, "Q");

    graphRealData->Draw("L same");

    TLegend *legend = new TLegend(0.1, 0.8, 0.3, 0.9);
    legend->AddEntry(graphRealData, "Real", "l");
    legend->AddEntry(f1, "Reconstruct", "l");

    legend->Draw();

#if configEnableDebug
    constexpr double conversionFactor = TMath::Ccgs() * 1e-9; // conversion factor from ns/cm to 1/c
    const double betaReciprocal = f1->GetParameter(0) * conversionFactor;
    printf("[Info] The reconstructed 1 / beta using the linear method: %f\n", betaReciprocal);
    printf("[Info] The real 1 / beta: %f\n", 1 / particle.getBeta());
#endif

    Detector_plotReconstructDataCanvas->SaveAs(fileName.c_str());
}

void Detector::plotDistributionOfReconstructionUsingLinearMethod(const Particle &particle, const int nReconstructions, const std::string &fileName) const
{
    const std::vector<std::tuple<double, double, TVector3>> hitData = this->particleHitData(particle, true);
    const int n = hitData.size();
    std::vector<double> hitTimes, propagationLengths;
    hitTimes.reserve(n), propagationLengths.reserve(n);

    for (const auto &hit : hitData)
    {
        hitTimes.push_back(std::get<0>(hit));
        propagationLengths.push_back(std::get<1>(hit));
    }

    const double betaReciprocalReal = 1 / particle.getBeta();

    TH1F *histogram = new TH1F("1/#beta_{real} - 1/#beta_{rec}", ";1/#beta_{real} - 1/#beta_{rec};Counts", 100, -0.1, 0.1);
    for (int i = 0; i < nReconstructions; ++i)
    {
        const std::vector<double> detectedTimes = this->detect(hitTimes);
        const double betaReciprocalReconstruction = this->reconstructUsingLinearMethod(detectedTimes, propagationLengths);
        histogram->Fill(betaReciprocalReal - betaReciprocalReconstruction);
    }

    TCanvas *Detector_plotDistributionOfReconstructionUsingLinearMethodCanvas = new TCanvas("Detector_plotDistributionOfReconstructionUsingLinearMethodCanvas", "", 3508, 2480); // A4 size in pixels(300 dpi)
    Detector_plotDistributionOfReconstructionUsingLinearMethodCanvas->cd();

    histogram->Draw();
    histogram->Fit("gaus");

    Detector_plotDistributionOfReconstructionUsingLinearMethodCanvas->SaveAs(fileName.c_str());
}