#include "Detector.h"

#include <algorithm>

#include <TMath.h>

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
#if configEnableDebugAtCoding
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
