#ifndef PARTICLE_H
#define PARTICLE_H

#include <Math/Point3D.h>
#include <Math/Vector3D.h>

#include "Config.h"

class Particle {
public:
  /* BEGIN Constructor & Destructor */

  /**
   * @brief Constructor
   *
   * The default constructor sets every member to zero, do NOT use it if you do not know what you are doing!
   */
  Particle();

  /**
   * @brief Constructor
   * @param charge Charge of the particle in e
   * @param mass0 Rest mass of the particle in MeV/c^2
   * @param momentum Momentum of the particle in MeV/c
   * @param position Position of the particle in cm
   */
  Particle(
      const int                    charge,
      const double                 mass0,
      const ROOT::Math::XYZVector &momentum,
      const ROOT::Math::XYZPoint  &position
  );

  /**
   * @brief Constructor
   * @param charge Charge of the particle in e
   * @param mass0 Rest mass of the particle in MeV/c^2
   * @param beta Beta of the particle (v/c)
   * @param position Position of the particle in cm
   * @param direction Moving direction of the particle
   */
  Particle(
      const int                    charge,
      const double                 mass0,
      double                       beta,
      const ROOT::Math::XYZPoint  &position,
      const ROOT::Math::XYZVector &direction = ROOT::Math::XYZVector(0, 0, 1)
  );

  ~Particle();

  /* END Constructor & Destructor */

  /* BEGIN Getters */

  /**
   * @brief Get the charge of the particle
   */
  inline int getCharge() const;

  /**
   * @brief Get the rest mass of the particle
   */
  inline double getMass0() const;

  /**
   * @brief Get the mass of the particle
   */
  inline double getMass() const;

  /**
   * @brief Get the energy of the particle
   */
  inline double getEnergy() const;

  /**
   * @brief Get the Lorentz factor of the particle
   */
  inline double getGamma() const;

  /**
   * @brief Get the beta of the particle
   */
  inline double getBeta() const;

  /**
   * @brief Get the position of the particle
   */
  ROOT::Math::XYZPoint getPosition() const;

  /**
   * @brief Get the momentum of the particle
   */
  ROOT::Math::XYZVector getMomentum() const;

  /**
   * @brief Get the velocity of the particle
   */
  ROOT::Math::XYZVector getVelocity() const;

  /* END Getters */

  /* BEGIN Setters */

  /**
   * @brief Set the mass of the particle
   * @param mass Mass of the particle in MeV/c^2
   * @return True if the mass is set successfully, false otherwise
   */
  bool setMass(double mass);

  /**
   * @brief Set the energy of the particle
   * @param energy Energy of the particle in MeV
   * @return True if the energy is set successfully, false otherwise
   */
  inline bool setEnergy(const double energy);

  /**
   * @brief Set the beta of the particle
   * @param beta Beta of the particle (v/c)
   * @return True if the beta is set successfully, false otherwise
   */
  bool setBeta(double beta);

  /**
   * @brief Set the beta * gamma of the particle
   * @param betaGamma Beta * gamma of the particle
   * @return True if the beta * gamma is set successfully, false otherwise
   */
  bool setBetaGamma(double betaGamma);

  /**
   * @brief Set the position of the particle
   * @param position Position of the particle in cm
   * @return True if the position is set successfully, false otherwise
   */
  bool setPosition(const ROOT::Math::XYZPoint &position);

  /**
   * @brief Set the momentum of the particle
   * @param momentum Momentum of the particle in MeV/c
   * @return True if the momentum is set successfully, false otherwise
   */
  bool setMomentum(const ROOT::Math::XYZVector &momentum);

  /**
   * @brief Set the velocity of the particle
   * @param velocity Velocity of the particle in c
   * @return True if the velocity is set successfully, false otherwise
   */
  bool setVelocity(const ROOT::Math::XYZVector &velocity);

  /**
   * @brief Set the moving direction of the particle
   * @param direction Moving direction of the particle
   * @return True if the moving direction is set successfully, false otherwise
   */
  bool setDirection(ROOT::Math::XYZVector direction);

  /* END Setters */

  /* BEGIN Methods */

  /**
   * @brief Calculate the maximum possible energy transfer to an electron in a single collision (in MeV)
   *
   * NOTE! This can be used ONLY for a point-like particle with rest mass much greater than the electron mass!
   * @return Maximum possible energy transfer to an electron in a single collision (in MeV)
   */
  double Wmax() const;

  /* END Methods */

protected:
  int                   pCharge;    // Charge of the particle in e
  double                pMass0;     // Rest mass of the particle in MeV/c^2
  double                pMass;      // Mass of the particle in MeV/c^2
  double                pEnergy;    // Energy of the particle in MeV (energy = mass in natural units)
  double                pGamma;     // Lorentz factor of the particle
  double                pBeta;      // Beta of the particle
  ROOT::Math::XYZVector pMomentum;  // Momentum of the particle in MeV/c
  ROOT::Math::XYZVector pVelocity;  // Velocity of the particle in c
  ROOT::Math::XYZPoint  pPosition;  // Position of the particle in cm
  ROOT::Math::XYZVector pDirection; // Moving direction of the particle
};

inline int Particle::getCharge() const { return this->pCharge; }

inline double Particle::getMass0() const { return this->pMass0; }

inline double Particle::getMass() const { return this->pMass; }

inline double Particle::getEnergy() const { return this->pEnergy; }

inline double Particle::getGamma() const { return this->pGamma; }

inline double Particle::getBeta() const { return this->pBeta; }

inline ROOT::Math::XYZPoint Particle::getPosition() const { return this->pPosition; }

inline ROOT::Math::XYZVector Particle::getMomentum() const { return this->pMomentum; }

inline ROOT::Math::XYZVector Particle::getVelocity() const { return this->pVelocity; }

inline bool Particle::setEnergy(const double energy) { return this->setMass(energy); }

#endif /* PARTICLE_H */