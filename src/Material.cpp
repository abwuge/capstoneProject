#include "Material.h"

#include <TMath.h>

Material::Material(
    const double density,
    const double Z,
    const double A,
    const double I,
    const double hbarOmegaP,
    const double DECa,
    const double DECk,
    const double DECx0,
    const double DECx1,
    const double DECdelta0,
    const double DECoverlineC
)
    : density(density), Z(Z), A(A), I(I), hbarOmegaP(hbarOmegaP), DECa(DECa), DECk(DECk), DECx0(DECx0), DECx1(DECx1),
      DECdelta0(DECdelta0), DECoverlineC(DECoverlineC) {}

Material::~Material() {}

double Material::delta(const double beta, const double gamma) const {
  const double x = TMath::Log10(beta * gamma);
  if (x <= this->DECx0) {
    return this->DECdelta0 ? this->DECdelta0 * TMath::Power(10, 2 * (x - this->DECx0)) : 0;
  } else if (x <= this->DECx1) {
    return 2 * TMath::Log(10) * x - this->DECoverlineC + this->DECa * TMath::Power(this->DECx1 - x, this->DECk);
  } else {
    return 2 * TMath::Log(10) * x - this->DECoverlineC;
  }
}

double Material::massStoppingPower(const Particle &particle) const {
  constexpr double K            = 0.307075; // MeV mol^-1 cm^2 (4 * pi * N_A * r_e^2 * m_e * c^2, coefficient for dE/dx)
  constexpr double electronMass = 0.51099895000; // Rest mass of the electron in MeV/c^2
  const int        z            = particle.getCharge(
  ); // charge number of the particle (since we use charge in units of e, the charge number is the same as the charge)
  const double beta  = particle.getBeta();
  const double gamma = particle.getGamma();

  if (Config::enableWarning) {
    const double betaGamma = beta * gamma;
    if (betaGamma < 0.1)
      printf(
          "[Warning] The energy of the particle (%g MeV/c) is too low! The error of the mass stopping power may be "
          "high!\n",
          particle.getEnergy()
      );
    else if (betaGamma > 1000)
      printf(
          "[Warning] The energy of the particle (%g MeV/c) is too high! The error of the mass stopping power may be "
          "high!\n",
          particle.getEnergy()
      );
  }

  const double beta2 = beta * beta;

  // <-dE/dx> = partA * [0.5 * ln(partB) - beta2 - 0.5 * delta(beta, gamma)]
  const double partA = K * (z * z) * (this->Z / this->A) * (1 / beta2);
  const double partB = 1e12 * (2 * electronMass * beta * beta * gamma * gamma * particle.Wmax())
                     / (this->I * this->I); // electronMass in MeV/c^2, Wmax in MeV, but, I in eV, so multiply by 1e12

  if (Config::enableDebug) printf("[Info] W_max: %g, partA: %g, partB: %g\n", particle.Wmax(), partA, partB);

  return partA * (0.5 * TMath::Log(partB) - beta2 - 0.5 * this->delta(beta, gamma));
}

double Material::linearStoppingPower(const Particle &particle) const {
  return this->density * this->massStoppingPower(particle);
}
