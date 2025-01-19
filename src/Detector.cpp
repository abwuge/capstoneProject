#include "Detector.h"

#include <algorithm>
#include <thread>

#include <TMath.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TString.h>
#include <TF1.h>
#include <TLine.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>

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
#if configUseBetheBloch
                const double particleEnergyLoss = scintillatorCounter.energyLoss(particle);
#else
                const double particleEnergyLoss = scintillatorCounter.LandauMostProbableEnergyLoss(particle);
#endif
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
    Detector_plotDeltaTimeCanvas->SetGrid();
    Detector_plotDeltaTimeCanvas->cd();

    TGraph *graph = new TGraph(this->scintillatorCounters.size());
    graph->SetTitle(";L [cm];#Deltat [ns]");

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

    // Zoom in (Use "{}" to limit the scope of variables)
    {
        double *x = graph->GetX(), *y = graph->GetY();
        int n = graph->GetN();
        double dx = x[n - 1] - x[0], dy = y[n - 1] - y[0];

        if (n > 1)
        {
            graph->GetXaxis()->SetRangeUser(x[1] - 1e-2 * dx, x[1] + 1e-2 * dx);
            graph->GetYaxis()->SetRangeUser(y[1] - 1e-2 * dy, y[1] + 1e-2 * dy);
            Detector_plotDeltaTimeCanvas->SaveAs(TString::Format("%s_zoom_1.png", fileName.c_str()));
        }

        if (n > 2)
        {
            graph->GetXaxis()->SetRangeUser(x[2] - 1e-2 * dx, x[2] + 1e-2 * dx);
            graph->GetYaxis()->SetRangeUser(y[2] - 1e-2 * dy, y[2] + 1e-2 * dy);
            Detector_plotDeltaTimeCanvas->SaveAs(TString::Format("%s_zoom_2.png", fileName.c_str()));
        }
    }
}

std::vector<double> Detector::detect(const Particle &particle) const
{
    std::vector<std::tuple<double, double, TVector3>> hitData = this->particleHitData(particle, true);
    std::vector<double> detectedTimes;
    detectedTimes.reserve(hitData.size());

    for (int i = 0; i < hitData.size(); ++i)
    {
        const double hitTime = std::get<0>(hitData.at(i));
#if configEnableTimeResolution
        detectedTimes.push_back(this->Random->Gaus(hitTime, this->scintillatorCounters.at(i).getTimeResolution()));
#else
        detectedTimes.push_back(hitTime);
#endif
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
#if configEnableTimeResolution
        detectedTimes.push_back(this->Random->Gaus(hitTime, this->scintillatorCounters.at(i).getTimeResolution()));
#else
        detectedTimes.push_back(hitTime);
#endif
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

    delete graph;
    delete f1;

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
#if configEnableDebug
        printf("[Info] propagation length: %f cm, hit time: %f ns\n", std::get<1>(hit), std::get<0>(hit));
#endif
    }

    std::vector<double> detectedTimes = this->detect(hitTimes);

    TGraph *graphRealData = new TGraph(n, &propagationLengths[0], &hitTimes[0]);
    graphRealData->SetMarkerStyle(kFullTriangleUp);
    graphRealData->SetMarkerSize(3);
    graphRealData->SetTitle(";L [cm];t [ns]");

    TGraph *graphReconstructData = new TGraph(n, &propagationLengths[0], &detectedTimes[0]);
    graphReconstructData->SetMarkerStyle(kFullCircle);
    graphReconstructData->SetMarkerColor(kRed);
    graphReconstructData->SetMarkerSize(2);

    TCanvas *Detector_plotReconstructDataCanvas = new TCanvas("Detector_plotReconstructDataCanvas", "", 3508, 2480); // A4 size in pixels(300 dpi)
    Detector_plotReconstructDataCanvas->SetGrid();
    Detector_plotReconstructDataCanvas->cd();

    TF1 *f1 = new TF1("f1", "[0] * x", 0, propagationLengths.back());

    graphRealData->Draw("AP");

    graphReconstructData->Draw("PSAME");
    graphReconstructData->Fit(f1, "Q");

    TLegend *legend = new TLegend(0.2, 0.75, 0.5, 0.85);
    legend->SetBorderSize(kNone);
    legend->AddEntry(graphRealData, "Real", "p");
    legend->AddEntry(graphReconstructData, "Detected", "p");
    const int exponent = TMath::FloorNint(TMath::Log10(f1->GetParameter(0)));
    const double mantissa = f1->GetParameter(0) / TMath::Power(10, exponent);
    legend->AddEntry(f1, Form("y = %.3f#times10^{%i} x (fit line)", mantissa, exponent));

    legend->Draw();

#if configEnableDebug
    constexpr double conversionFactor = TMath::Ccgs() * 1e-9; // conversion factor from ns/cm to 1/c
    const double betaReciprocal = f1->GetParameter(0) * conversionFactor;
    printf("[Info] The reconstructed 1 / beta using the linear method: %f\n", betaReciprocal);
    printf("[Info] The real 1 / beta: %f\n", 1 / particle.getBeta());
#endif

    Detector_plotReconstructDataCanvas->SaveAs(fileName.c_str());
}

void Detector::processReconstruction(const Particle &particle, const bool enableLinerMethod, const double betaReciprocalReal, TH1F *histogram) const
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

    double betaReciprocalReconstruction;
    if (enableLinerMethod)
        betaReciprocalReconstruction = this->reconstructUsingLinearMethod(this->detect(hitTimes), propagationLengths);
    else
        betaReciprocalReconstruction = this->reconstructUsingNonLinearMethod(particle, this->detect(hitTimes), propagationLengths);

    histogram->Fill(betaReciprocalReal - betaReciprocalReconstruction);
}

std::pair<double, double> Detector::distributionOfReconstruction(const Particle &particle, const int nReconstructions, const bool enableLinerMethod, const bool enablePlot, const std::string &fileName) const
{
    const double betaReciprocalReal = 1 / particle.getBeta();

    TH1F *histogram = new TH1F("#Delta(1/#beta)", ";#Delta(1/#beta);", 100, -0.1, 0.1);

    std::vector<std::future<void>> futures;
    futures.reserve(nReconstructions);
    for (int i = 0; i < nReconstructions; ++i)
        futures.push_back(ThreadPool::getInstance().enqueue([this, &particle, enableLinerMethod, betaReciprocalReal, histogram]()
                                                            { this->processReconstruction(particle, enableLinerMethod, betaReciprocalReal, histogram); }));

    for (auto &future : futures)
        future.get();

    TCanvas *Detector_plotDistributionOfReconstructionCanvas;
    if (enablePlot)
    {
        Detector_plotDistributionOfReconstructionCanvas = new TCanvas("Detector_plotDistributionOfReconstructionCanvas", "", 3508, 2480); // A4 size in pixels(300 dpi)
        Detector_plotDistributionOfReconstructionCanvas->cd();

        histogram->Draw();

        TFitResultPtr fitResult = histogram->Fit("gaus", "QS");

        double mean = fitResult->Value(1), sigma = fitResult->Value(2);

        // TLine *line = new TLine(deltaBetaReciprocalWithZeroResolution, 0, deltaBetaReciprocalWithZeroResolution, histogram->GetMaximum());
        // line->Draw("same");

        TLegend *legend = new TLegend(0.6, 0.75, 0.85, 0.85);
        legend->SetBorderSize(kNone);
        legend->SetHeader(Form("#mu = %.3f, #sigma = %.3f", mean, sigma), "C");
        // legend->AddEntry(line, Form("#Delta(1/#beta)_{real} = %.3f", deltaBetaReciprocalWithZeroResolution), "l");
        legend->Draw();

#if configEnableDebug
        printf("[Info] The mean of the distribution of the difference between real and reconstructed 1/beta: %f\n", mean);
        // printf("[Info] The difference between real and reconstructed 1/beta with zero resolution: %f\n", deltaBetaReciprocalWithZeroResolution);
#endif

        Detector_plotDistributionOfReconstructionCanvas->SaveAs(fileName.c_str());
    }

    std::pair<double, double> meanAndStandardDeviation = {histogram->GetMean(), histogram->GetRMS()};

    delete histogram;
    delete Detector_plotDistributionOfReconstructionCanvas;

    return meanAndStandardDeviation;
}

void Detector::plotDeltaBetaReciprocal(Particle particle, const double betaMin, const double betaMax, const int nPoints, const bool enableLinearMethod, const std::string &fileName) const
{
    TCanvas *Detector_plotDeltaBetaReciprocalCanvas = new TCanvas("Detector_plotDeltaBetaReciprocalCanvas", "", 3508, 2480); // A4 size in pixels(300 dpi)
    Detector_plotDeltaBetaReciprocalCanvas->SetGrid();
    Detector_plotDeltaBetaReciprocalCanvas->cd();

    TGraphErrors *graphErrors = new TGraphErrors(nPoints + 1);
    graphErrors->SetTitle(";#beta;#Delta(1/#beta)");

    const double step = (betaMax - betaMin) / nPoints;
    for (int i = 0; i <= nPoints; ++i)
    {
        const double beta = betaMin + i * step;
        particle.setBeta(beta);
        const std::pair<double, double> deltaBetaReciprocal = this->distributionOfReconstruction(particle, 10000, enableLinearMethod, false);
        graphErrors->SetPoint(i, beta, deltaBetaReciprocal.first);
        graphErrors->SetPointError(i, 0, deltaBetaReciprocal.second);
    }

    graphErrors->Draw("AL");

    Detector_plotDeltaBetaReciprocalCanvas->SaveAs(fileName.c_str());
}

double Detector::reconstructUsingNonLinearMethod(const Particle &particle) const
{
    const std::vector<std::tuple<double, double, TVector3>> hitData = this->particleHitData(particle, true);
    std::vector<double> hitTimes, propagationLengths;
    hitTimes.reserve(hitData.size()), propagationLengths.reserve(hitData.size());

    for (const auto &hit : hitData)
    {
        hitTimes.push_back(std::get<0>(hit));
        propagationLengths.push_back(std::get<1>(hit));
    }

    return this->reconstructUsingNonLinearMethod(particle, this->detect(hitTimes), propagationLengths);
}

double Detector::reconstructUsingNonLinearMethod(Particle particle, const std::vector<double> &detectedTimes, const std::vector<double> &propagationLengths) const
{
    const double initialBetaReciprocal = this->reconstructUsingLinearMethod(detectedTimes, propagationLengths);

    auto chi2Function = [&](const double *parameters)
    {
        particle.setBeta(1 / parameters[0]);

        std::vector<std::tuple<double, double, TVector3>> reconstructedHitData = this->particleHitData(particle, true);
        std::vector<double> reconstructedHitTimes;
        reconstructedHitTimes.reserve(reconstructedHitData.size());

        for (const auto &hit : reconstructedHitData)
            reconstructedHitTimes.push_back(std::get<0>(hit));

        double chi2 = 0;
        for (size_t i = 0; i < detectedTimes.size(); ++i)
        {
            double timeDiff = reconstructedHitTimes.at(i) - detectedTimes.at(i);
            double timeResolution = this->getScintillatorCounters().at(i).getTimeResolution();
            chi2 += TMath::Power(timeDiff, 2) / TMath::Power(timeResolution, 2);
        }
        return chi2;
    };

    ROOT::Math::Minimizer *minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

    ROOT::Math::Functor functor(chi2Function, 1);
    minimizer->SetFunction(functor);

    minimizer->SetLimitedVariable(0, "betaReciprocal", initialBetaReciprocal, 1e-5, 1, 10);

#if configEnableDebug
    minimizer->SetPrintLevel(1);
#endif

    minimizer->Minimize();
    const double betaReciprocal = minimizer->X()[0];

#if configEnableDebug
    printf("[Info] The reconstructed 1/beta using the linear method: %f\n", initialBetaReciprocal);
    printf("[Info] The reconstructed 1/beta using the non-linear method: %f\n", betaReciprocal);
#endif

    return betaReciprocal;
}
