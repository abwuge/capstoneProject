#include "Particle.h"

#include <TMath.h>

Particle::Particle()
    : charge(0), mass0(0), mass(0), energy(0), gamma(0), beta(0), position(ROOT::Math::XYZVector(0, 0, 0)),
      momentum(ROOT::Math::XYZVector(0, 0, 0)), velocity(ROOT::Math::XYZVector(0, 0, 0)) {}

Particle::Particle(
    const int                    charge,
    const double                 mass0,
    const ROOT::Math::XYZVector &momentum,
    const ROOT::Math::XYZVector &position
)
    : charge(charge), mass0(mass0), momentum(momentum), position(position) {
  this->energy   = TMath::Sqrt(this->momentum.Mag2() + this->mass0 * this->mass0);
  this->mass     = this->energy;
  this->gamma    = this->energy / this->mass0;
  this->velocity = this->momentum * (1 / this->mass);
  this->beta     = this->velocity.R();
}

Particle::Particle(
    const int                    charge,
    const double                 mass0,
    double                       beta,
    const ROOT::Math::XYZVector &position,
    const ROOT::Math::XYZVector &direction
)
    : charge(charge), mass0(mass0), position(position) {
  if (beta < 0) {
    if (Config::enableWarning) printf("[Warning] Beta cannot be negative, setting beta to 0!\n");

    beta = 0;
  }

  if (beta >= 1) {
    if (Config::enableWarning) printf("[Warning] Beta must be less than 1, setting beta to 0.9999999999999999!\n");

    beta = 0.9999999999999999;
  }

  this->beta     = beta;
  this->gamma    = 1 / TMath::Sqrt(1 - this->beta * this->beta);
  this->energy   = this->gamma * this->mass0;
  this->mass     = this->energy;
  this->velocity = this->beta * direction;
  this->momentum = this->mass * this->velocity;
}

Particle::~Particle() {}

bool Particle::setMass(double mass) {
  if (this->momentum.R() == 0) {
    if (Config::enableWarning)
      printf("[Warning] No direction is set for the particle! Use (0, 0, 1) as the direction!\n");

    this->momentum = ROOT::Math::XYZVector(0, 0, 1);
  }

  if (mass < this->mass0) {
    if (Config::enableWarning)
      printf("[Warning] Mass / Energy cannot be less than the rest mass! Mass / Energy is set to the rest mass!\n");

    mass = this->mass0;
  }

  this->mass     = mass;
  this->energy   = this->mass;
  this->gamma    = this->energy / this->mass0;
  this->momentum = TMath::Sqrt(this->energy * this->energy - this->mass0 * this->mass0) * this->momentum.Unit();
  this->velocity = this->momentum * (1 / this->mass);
  this->beta     = this->velocity.R();

  return true;
}

bool Particle::setBeta(double beta) {
  if (this->velocity.R() == 0) {
    if (Config::enableWarning)
      printf("[Warning] No direction is set for the particle! Use (0, 0, 1) as the direction!\n");

    this->velocity = ROOT::Math::XYZVector(0, 0, 1);
  }

  if (beta < 0 || beta > 1) {
    if (Config::enableWarning) printf("[Warning] Beta must be between 0 and 1! Beta remains %f!\n", this->beta);

    return false;
  }

  if (beta == 1) {
    if (Config::enableWarning) printf("[Warning] Beta cannot be 1! Beta is set to 0.9999999999999999!\n");

    beta = 0.9999999999999999;
  }

  this->beta     = beta;
  this->gamma    = 1 / TMath::Sqrt(1 - this->beta * this->beta);
  this->energy   = this->gamma * this->mass0;
  this->mass     = this->energy;
  this->velocity = this->beta * this->velocity.Unit();
  this->momentum = this->mass * this->velocity;

  return true;
}

bool Particle::setBetaGamma(double betaGamma) {
  if (this->velocity.R() == 0) {
    if (Config::enableWarning)
      printf("[Warning] No direction is set for the particle! Use (0, 0, 1) as the direction!\n");

    this->momentum = ROOT::Math::XYZVector(0, 0, 1);
  }

  if (betaGamma < 0) {
    if (Config::enableWarning) printf("[Warning] Beta * gamma cannot be negative! Beta * gamma is set to 0!\n");

    betaGamma = 0;
  }

  this->beta     = betaGamma / TMath::Sqrt(1 + betaGamma * betaGamma);
  this->gamma    = betaGamma / this->beta;
  this->energy   = this->gamma * this->mass0;
  this->mass     = this->energy;
  this->velocity = this->beta * this->velocity.Unit();
  this->momentum = this->mass * this->velocity;

  return true;
}

bool Particle::setPosition(const ROOT::Math::XYZVector &position) {
  this->position = position;

  return true;
}

bool Particle::setMomentum(const ROOT::Math::XYZVector &momentum) {
  this->momentum = momentum;
  this->energy   = TMath::Sqrt(this->momentum.Mag2() + this->mass0 * this->mass0);
  this->mass     = this->energy;
  this->gamma    = this->energy / this->mass0;
  this->velocity = this->momentum * (1 / this->mass);
  this->beta     = this->velocity.R();

  return true;
}

bool Particle::setVelocity(const ROOT::Math::XYZVector &velocity) {
  this->velocity = velocity;
  this->beta     = this->velocity.R();
  this->gamma    = 1 / TMath::Sqrt(1 - this->beta * this->beta);
  this->energy   = this->gamma * this->mass0;
  this->mass     = this->energy;
  this->momentum = this->mass * this->velocity;

  return true;
}

double Particle::Wmax() const {
  constexpr double electronMass = 0.51099895000; // Rest mass of the electron in MeV/c^2

  if (Config::enableWarning)
    if (this->mass0 < 100 * electronMass)
      printf(
          "[Warning] The mass of the particle is NOT santisfied the condition: M >> m_e! The error of W_max may be "
          "high!\n"
      );

  const double betaGamma = this->beta * this->gamma; // Beta * gamma of the particle
  const double electronMassMassRatio =
      electronMass / this->mass0;                    // Rest mass of the electron / rest mass of the particle

  return 2 * electronMass * betaGamma * betaGamma
       / (1 + 2 * this->gamma * electronMassMassRatio + electronMassMassRatio * electronMassMassRatio);
}
