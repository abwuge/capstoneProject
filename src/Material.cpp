#include "Material.h"

#include <TMath.h>

Material::Material(const double density, const double Z, const double A, const double I, const double hbarOmegaP,
                   const double DECa, const double DECk, const double DECx0, const double DECx1, const double DECdelta0, const double DECoverlineC)
    : density(density), Z(Z), A(A), I(I), hbarOmegaP(hbarOmegaP), DECa(DECa), DECk(DECk), DECx0(DECx0), DECx1(DECx1), DECdelta0(DECdelta0), DECoverlineC(DECoverlineC) {}

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
    const double x = TMath::Log10(beta * gamma);
    if (x <= DECx0)
    {
        return DECdelta0 ? DECdelta0 * TMath::Power(10, 2 * (x - DECx0)) : 0;
    }
    else if (x <= DECx1)
    {
        return 2 * TMath::Log(10) * x - DECoverlineC + DECa * TMath::Power(DECx1 - x, DECk);
    }
    else
    {
        return 2 * TMath::Log(10) * x - DECoverlineC;
    }
}

double Material::massStoppingPower(const Particle &particle) const
{
    constexpr double K = 0.307075;                 // MeV mol^-1 cm^2 (4 * pi * N_A * r_e^2 * m_e * c^2, coefficient for the Bethe formula)
    constexpr double electronMass = 0.51099895000; // Rest mass of the electron in MeV/c^2
    const double z = particle.getCharge();
    const double beta = particle.getBeta();
    const double gamma = particle.getGamma();

    const double beta2 = beta * beta;

    // <-dE/dx> = partA * [0.5 * ln(partB) - beta2 - 0.5 * delta(beta, gamma)]
    const double partA = K * (z * z) * (Z / A) * (1 / beta2);
    const double partB = 1e12 * (2 * electronMass * beta * beta * gamma * gamma * particle.Wmax()) / (I * I); // electronMass in MeV/c^2, Wmax in MeV, but, I in eV, so multiply by 1e12

    #ifdef DEBUG
        printf("W_max: %g, partB: %g\n", particle.Wmax(), partB);
    #endif

    return partA * (0.5 * TMath::Log(partB) - beta2 - 0.5 * delta(beta, gamma));
}

double Material::linearStoppingPower(const Particle &particle) const
{
    return density * massStoppingPower(particle);
}
