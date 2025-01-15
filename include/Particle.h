#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>

#include <TVector3.h>

#include "config.h"

class Particle
{
private:
    double charge; // Charge of the particle in e
    double mass0;  // Rest mass of the particle in MeV/c^2
    double mass;   // Mass of the particle in MeV/c^2
    // Question: I have used `double &energy = mass` to ensure that the energy is always equal to the mass, but it is not working as expected. Why?
    double energy; // Energy of the particle in MeV (energy = mass in natural units)
    double gamma;  // Lorentz factor of the particle
    double beta;   // Beta of the particle
    TVector3 position; // Position of the particle in cm
    TVector3 momentum; // Momentum of the particle in MeV/c
    TVector3 velocity; // Velocity of the particle in c

public:
    /* BEGIN Constructor & Destructor */

    /**
     * @brief Constructor
     * @param charge Charge of the particle in e
     * @param mass0 Rest mass of the particle in MeV/c^2
     * @param momentum Momentum of the particle in MeV/c
     * @param position Position of the particle in cm
     */
    Particle(const double charge, const double mass0, const TVector3 &momentum, const TVector3 &position);

    /**
     * @brief Constructor
     * @param charge Charge of the particle in e
     * @param mass0 Rest mass of the particle in MeV/c^2
     * @param beta Beta of the particle (v/c)
     * @param position Position of the particle in cm
     * @param direction Moving direction of the particle
     */
    Particle(const double charge, const double mass0, double beta, const TVector3 &position, const TVector3 &direction = TVector3(0, 0, 1));

    ~Particle();

    /* END Constructor & Destructor */

    /* BEGIN Getters */

    /**
     * @brief Get the charge of the particle
     */
    inline double getCharge() const;

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
    TVector3 getPosition() const;

    /**
     * @brief Get the momentum of the particle
     */
    TVector3 getMomentum() const;

    /**
     * @brief Get the velocity of the particle
     */
    TVector3 getVelocity() const;

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
    bool setBeta(const double beta);

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
    bool setPosition(const TVector3 &position);

    /**
     * @brief Set the momentum of the particle
     * @param momentum Momentum of the particle in MeV/c
     * @return True if the momentum is set successfully, false otherwise
     */
    bool setMomentum(const TVector3 &momentum);

    /**
     * @brief Set the velocity of the particle
     * @param velocity Velocity of the particle in c
     * @return True if the velocity is set successfully, false otherwise
     */
    bool setVelocity(const TVector3 &velocity);

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
};

inline double Particle::getCharge() const { return this->charge; }

inline double Particle::getMass0() const { return this->mass0; }

inline double Particle::getMass() const { return this->mass; }

inline double Particle::getEnergy() const { return this->energy; }

inline double Particle::getGamma() const { return this->gamma; }

inline double Particle::getBeta() const { return this->beta; }

inline TVector3 Particle::getPosition() const { return this->position; }

inline TVector3 Particle::getMomentum() const { return this->momentum; }

inline TVector3 Particle::getVelocity() const { return this->velocity; }

inline bool Particle::setEnergy(const double energy) { return this->setMass(energy); }

#endif /* PARTICLE_H */
