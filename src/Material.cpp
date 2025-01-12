#include "Material.h"
#include <TMath.h>

Material::Material(const double density, const double Z, const double A, const double zToARatio, const double I, const double hbarOmegaP)
    : density(density), Z(Z), A(A), zToARatio(zToARatio), I(I), hbarOmegaP(hbarOmegaP) {}

Material::~Material() {}

double Material::getDensity() const
{
    return density;
}

double Material::getZ() const
{
    return Z;
}

double Material::getA() const
{
    return A;
}

double Material::getZToARatio() const
{
    return zToARatio;
}

double Material::getI() const
{
    return I;
}

double Material::getHbarOmegaP() const
{
    return hbarOmegaP;
}

double Material::delta(const double beta, const double gamma) const
{
    return 2 * (TMath::Log(hbarOmegaP / I) - TMath::Log(beta * gamma) - 0.5);
}

double Material::meanRateOfEnergyLoss(const Particle &particle) const
{
    constexpr double K = 0.307075;                 // MeV mol^-1 cm^2 (4 * pi * N_A * r_e^2 * m_e * c^2, coefficient for the Bethe formula)
    constexpr double electronMass = 0.51099895000; // Rest mass of the electron in MeV/c^2
    const double z = particle.getCharge();
    const double beta = particle.getBeta();
    const double gamma = particle.getGamma();

    return density *
           (K * (z * z) * (Z / A) * (1 / (beta * beta)) *
            (0.5 * TMath::Log(2 * electronMass * beta * beta * gamma * gamma * particle.Wmax() / (I * I)) - beta * beta - 0.5 * delta(beta, gamma)));
}