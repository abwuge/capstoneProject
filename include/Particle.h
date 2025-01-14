#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <TVector3.h>

class Particle
{
private:
    double charge;     // Charge of the particle in e
    double mass0;      // Rest mass of the particle in MeV/c^2
    double mass;       // Mass of the particle in MeV/c^2
    double energy;     // Energy of the particle in MeV
    double gamma;      // Lorentz factor of the particle
    double beta;       // Beta of the particle
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
    double getCharge() const;

    /**
     * @brief Get the rest mass of the particle
     */
    double getMass0() const;

    /**
     * @brief Get the mass of the particle
     */
    double getMass() const;

    /**
     * @brief Get the energy of the particle
     */
    double getEnergy() const;

    /**
     * @brief Get the Lorentz factor of the particle
     */
    double getGamma() const;

    /**
     * @brief Get the beta of the particle
     */
    double getBeta() const;

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
     * @brief Set the beta of the particle
     * @param beta Beta of the particle (v/c)
     * @return True if the beta is set successfully, false otherwise
     */
    bool setBeta(const double beta);

    /**
     * @brief Set the momentum of the particle
     * @param momentum Momentum of the particle in MeV/c
     * @return True if the momentum is set successfully, false otherwise
     */
    bool setMomentum(const TVector3 &momentum);

    /**
     * @brief Set the beta * gamma of the particle
     * @param betaGamma Beta * gamma of the particle
     * @return True if the beta * gamma is set successfully, false otherwise
     */
    bool setBetaGamma(const double betaGamma);

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

#endif /* PARTICLE_H */
