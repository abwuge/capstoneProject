#include "Particle.h"

#include <TMath.h>

Particle::Particle()
    : pCharge(0), pMass0(0), pMass(0), pEnergy(0), pGamma(0), pBeta(0), pMomentum(ROOT::Math::XYZVector(0, 0, 0)),
      pVelocity(ROOT::Math::XYZVector(0, 0, 0)), pPosition(ROOT::Math::XYZPoint(0, 0, 0)),
      pDirection(ROOT::Math::XYZVector(0, 0, 1)) {}

Particle::Particle(
    const int                    charge,
    const double                 mass0,
    const ROOT::Math::XYZVector &momentum,
    const ROOT::Math::XYZPoint  &position
)
    : pCharge(charge), pMass0(mass0), pMomentum(momentum), pPosition(position) {
  this->pEnergy    = TMath::Sqrt(this->pMomentum.Mag2() + this->pMass0 * this->pMass0);
  this->pMass      = this->pEnergy;
  this->pGamma     = this->pEnergy / this->pMass0;
  this->pVelocity  = this->pMomentum * (1 / this->pMass);
  this->pBeta      = this->pVelocity.R();
  this->pDirection = this->pMomentum.Unit();
}

Particle::Particle(
    const int                    charge,
    const double                 mass0,
    double                       beta,
    const ROOT::Math::XYZPoint  &position,
    const ROOT::Math::XYZVector &direction
)
    : pCharge(charge), pMass0(mass0), pPosition(position), pDirection(direction.Unit()) {
  if (!this->pDirection.R()) {
    if (Config::enableWarning)
      printf("[Warning] The direction of the particle cannot be (0, 0, 0)! Setting the direction to (0, 0, 1)!\n");

    this->pDirection = ROOT::Math::XYZVector(0, 0, 1);
  }

  if (beta < 0) {
    if (Config::enableWarning) printf("[Warning] Beta cannot be negative, setting beta to 0!\n");

    beta = 0;
  }

  if (beta >= 1) {
    if (Config::enableWarning) printf("[Warning] Beta must be less than 1, setting beta to 0.9999999999999999!\n");

    beta = 0.9999999999999999;
  }

  this->pBeta     = beta;
  this->pGamma    = 1 / TMath::Sqrt(1 - this->pBeta * this->pBeta);
  this->pEnergy   = this->pGamma * this->pMass0;
  this->pMass     = this->pEnergy;
  this->pVelocity = this->pBeta * this->pDirection;
  this->pMomentum = this->pMass * this->pVelocity;
}

Particle::~Particle() {}

bool Particle::setMass(double mass) {
  if (mass < this->pMass0) {
    if (Config::enableWarning)
      printf("[Warning] Mass / Energy cannot be less than the rest mass! Mass / Energy is set to the rest mass!\n");

    mass = this->pMass0;
  }

  this->pMass     = mass;
  this->pEnergy   = this->pMass;
  this->pGamma    = this->pEnergy / this->pMass0;
  this->pMomentum = TMath::Sqrt(this->pEnergy * this->pEnergy - this->pMass0 * this->pMass0) * this->pDirection;
  this->pVelocity = this->pMomentum * (1 / this->pMass);
  this->pBeta     = this->pVelocity.R();

  return true;
}

bool Particle::setBeta(double beta) {
  if (beta < 0 || beta > 1) {
    if (Config::enableWarning) printf("[Warning] Beta must be between 0 and 1! Beta remains %f!\n", this->pBeta);

    return false;
  }

  if (beta > 0.9999999999999999) {
    if (Config::enableWarning) printf("[Warning] Beta cannot be 1! Beta is set to 0.9999999999999999!\n");

    beta = 0.9999999999999999;
  }

  this->pBeta     = beta;
  this->pGamma    = 1 / TMath::Sqrt(1 - this->pBeta * this->pBeta);
  this->pEnergy   = this->pGamma * this->pMass0;
  this->pMass     = this->pEnergy;
  this->pVelocity = this->pBeta * this->pDirection;
  this->pMomentum = this->pMass * this->pVelocity;

  return true;
}

bool Particle::setBetaGamma(double betaGamma) {
  if (betaGamma < 0) {
    if (Config::enableWarning) printf("[Warning] Beta * gamma cannot be negative! Beta * gamma is set to 0!\n");

    betaGamma = 0;
  }

  this->pBeta     = betaGamma / TMath::Sqrt(1 + betaGamma * betaGamma);
  this->pGamma    = betaGamma / this->pBeta;
  this->pEnergy   = this->pGamma * this->pMass0;
  this->pMass     = this->pEnergy;
  this->pVelocity = this->pBeta * this->pDirection;
  this->pMomentum = this->pMass * this->pVelocity;

  return true;
}

bool Particle::setPosition(const ROOT::Math::XYZPoint &position) {
  this->pPosition = position;

  return true;
}

bool Particle::setMomentum(const ROOT::Math::XYZVector &momentum) {
  this->pMomentum  = momentum;
  this->pDirection = this->pMomentum.Unit();
  this->pEnergy    = TMath::Sqrt(this->pMomentum.Mag2() + this->pMass0 * this->pMass0);
  this->pMass      = this->pEnergy;
  this->pGamma     = this->pEnergy / this->pMass0;
  this->pVelocity  = this->pMomentum * (1 / this->pMass);
  this->pBeta      = this->pVelocity.R();

  return true;
}

bool Particle::setVelocity(const ROOT::Math::XYZVector &velocity) {
  const double beta = velocity.R();
  if (beta >= 1) {
    printf("[Error] The velocity of the particle must be less than c!\n");

    return false;
  }

  this->pVelocity  = velocity;
  this->pBeta      = beta;
  this->pGamma     = 1 / TMath::Sqrt(1 - this->pBeta * this->pBeta);
  this->pEnergy    = this->pGamma * this->pMass0;
  this->pMass      = this->pEnergy;
  this->pMomentum  = this->pMass * this->pVelocity;
  this->pDirection = this->pMomentum.Unit();

  return true;
}

bool Particle::setDirection(ROOT::Math::XYZVector direction) {
  this->pDirection = direction.Unit();

  if (!this->pDirection.R()) {
    if (Config::enableWarning)
      printf("[Warning] The direction of the particle cannot be (0, 0, 0)! Setting the direction to (0, 0, 1)!\n");

    this->pDirection = ROOT::Math::XYZVector(0, 0, 1);
  }

  this->pMomentum = this->pMomentum.R() * this->pDirection;
  this->pVelocity = this->pVelocity.R() * this->pDirection;

  return true;
}

double Particle::Wmax() const {
  constexpr double electronMass = 0.51099895000; // Rest mass of the electron in MeV/c^2

  if (Config::enableWarning)
    if (this->pMass0 < 100 * electronMass)
      printf(
          "[Warning] The mass of the particle is NOT santisfied the condition: M >> m_e! The error of W_max may be "
          "high!\n"
      );

  const double betaGamma = this->pBeta * this->pGamma; // Beta * gamma of the particle
  const double electronMassMassRatio =
      electronMass / this->pMass0;                     // Rest mass of the electron / rest mass of the particle

  return 2 * electronMass * betaGamma * betaGamma
       / (1 + 2 * this->pGamma * electronMassMassRatio + electronMassMassRatio * electronMassMassRatio);
}