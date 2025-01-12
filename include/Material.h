#ifndef MATERIAL_H
#define MATERIAL_H

#include "Particle.h"

class Material
{
private:
    double density;    // Density of the material in g/cm^3
    double Z;          // Atomic number of the material
    double A;          // Atomic mass of the material in g/mol
    double zToARatio;  // <Z/A> of the material in g/mol
    double I;          // Mean excitation energy of the material in eV
    double hbarOmegaP; // Plasma energy of the material in eV

public:
    /* BEGIN Constructor & Destructor */

    /**
     * @brief Constructor
     * @param density Density of the material in g/cm^3
     * @param Z Atomic number of the material
     * @param A Atomic mass of the material in g/mol
     * @param zToARatio <Z/A> of the material in g/mol
     * @param I Mean excitation energy of the material in eV
     * @param hbarOmegaP Plasma energy of the material in eV
     */
    Material(const double density, const double Z, const double A, const double zToARatio, const double I, const double hbarOmegaP);

    ~Material();

    /* END Constructor & Destructor */

    /* BEGIN Getters */

    /**
     * @brief Get the density of the material in g/cm^3
     */
    double getDensity() const;

    /**
     * @brief Get the atomic number of the material
     */
    double getZ() const;

    /**
     * @brief Get the atomic mass of the material in g/mol
     */
    double getA() const;

    /**
     * @brief Get the <Z/A> of the material in g/mol
     */
    double getZToARatio() const;

    /**
     * @brief Get the mean excitation energy of the material in eV
     */
    double getI() const;

    /**
     * @brief Get the plasma energy of the material in eV
     */
    double getHbarOmegaP() const;

    /* END Getters */

    /* BEGIN Methods */

    /**
     * @brief Calculate the density effect correction to ionization energy loss
     *
     * NOTE! This can be used ONLY for energies at a very high level!
     * @param beta Velocity of the particle in units of c
     * @param gamma Lorentz factor of the particle
     * @return Density effect correction to ionization energy loss
     */
    double delta(const double beta, const double gamma) const;

    /**
     * @brief Calculate the mean rate of energy loss of the particle in the material
     * @param particle Incident particle
     * @return <-dE/dx>: Mean rate of energy loss of the particle in the material in MeV/cm
     */
    double meanRateOfEnergyLoss(const Particle &particle) const;

    /* END Methods */
};

#endif /* MATERIAL_H */