#ifndef MATERIAL_H
#define MATERIAL_H

#include <unordered_map>

#include "Config.h"
#include "Particle.h"

/**
 * @brief Enum class for material names
 *
 * The materials in the enum class is predefined materials
 */
enum class MaterialName
{
    Copper,
    Polyvinyltoluene,
    Polystyrene,
};

class Material
{
protected:
    double density;      // Density of the material in g/cm^3
    double Z;            // Atomic number of the material
    double A;            // Atomic mass of the material in g/mol
    double I;            // Mean excitation energy of the material in eV
    double hbarOmegaP;   // Plasma energy of the material in eV
    double DECa;         // Density effect correction parameter a
    double DECk;         // Density effect correction parameter k
    double DECx0;        // Density effect correction parameter \x_0
    double DECx1;        // Density effect correction parameter \x_1
    double DECdelta0;    // Density effect correction parameter \delta_0
    double DECoverlineC; // Density effect correction parameter \overline{C}

public:
    /* BEGIN Constructor & Destructor */

    /**
     * @brief Constructor
     * @param density Density of the material in g/cm^3
     * @param Z Atomic number of the material
     * @param A Atomic mass of the material in g/mol
     * @param I Mean excitation energy of the material in eV
     * @param hbarOmegaP Plasma energy of the material in eV
     * @param DECa Density effect correction parameter a
     * @param DECk Density effect correction parameter k
     * @param DECx0 Density effect correction parameter x_0
     * @param DECx1 Density effect correction parameter x_1
     * @param DECdelta0 Density effect correction parameter delta_0 (Put 0 if the material is nonconducting)
     * @param DECoverlineC Density effect correction parameter overline{C}
     */
    Material(const double density, const double Z, const double A, const double I, const double hbarOmegaP, const double DECa, const double DECk, const double DECx0, const double DECx1, const double DECdelta0, const double DECoverlineC);

    /**
     * @brief Copy constructor
     * @param material Material to copy
     */
    Material(const Material &material) = default;

    ~Material();

    /* END Constructor & Destructor */

    /* BEGIN Getters */

    /**
     * @brief Get the density of the material in g/cm^3
     */
    inline double getDensity() const;

    /**
     * @brief Get the atomic number of the material
     */
    inline double getZ() const;

    /**
     * @brief Get the atomic mass of the material in g/mol
     */
    inline double getA() const;

    /**
     * @brief Get the mean excitation energy of the material in eV
     */
    inline double getI() const;

    /**
     * @brief Get the plasma energy of the material in eV
     */
    inline double getHbarOmegaP() const;

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
     * @brief Calculate the mass stopping power of the particle in the material
     * @param particle Incident particle
     * @return <-dE/dx>: Mass stopping power of the particle in the material in MeV·cm^2/g
     */
    double massStoppingPower(const Particle &particle) const;

    /**
     * @brief Calculate the linear stopping power of the particle in the material
     * @param particle Incident particle
     * @return <-dE/dx>: Linear stopping power of the particle in the material in MeV/cm
     */
    double linearStoppingPower(const Particle &particle) const;

    /* END Methods */
};

inline double Material::getDensity() const { return this->density; }

inline double Material::getZ() const { return this->Z; }

inline double Material::getA() const { return this->A; }

inline double Material::getI() const { return this->I; }

inline double Material::getHbarOmegaP() const { return this->hbarOmegaP; }

const std::unordered_map<MaterialName, Material> materials = {
    /* TODO: IDK if delta(x_0) is delta_0 */
    {MaterialName::Copper, Material(8.960, 29, 63.546, 322.0, 58.27,              // from https://pdg.lbl.gov/2024/AtomicNuclearProperties/HTML/copper_Cu.html
                                    0.2557, 2.613, -0.0089, 3.0, 0.0893, 4.419)}, // from doi:10.1103/physrevb.26.6067
    /* TODO: I don't know the density effect correction parameters for Polyvinyltoluene */
    {MaterialName::Polyvinyltoluene, Material(1.032, 1 * 10.00 + 6 * 9.03, 1.0080 * 10.00 + 12.0107 * 9.03, 64.7, 21.54, // from https://pdg.lbl.gov/2024/AtomicNuclearProperties/HTML/polyvinyltoluene.html
                                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0)},                                            // I don't know the density effect correction parameters for Polyvinyltoluene
    /* This material is similar to polyvinyltoluene */
    {MaterialName::Polystyrene, Material(1.060, 1 * 8 + 6 * 8, 1.0080 * 8 + 12.0107 * 8, 68.7, 21.75, // from https://pdg.lbl.gov/2024/AtomicNuclearProperties/HTML/polystyrene.html
                                         0.3670, 2.724, 0.1647, 2.2, 0.0, 3.300)}                     // from doi:10.1103/physrevb.26.6067
};

#endif /* MATERIAL_H */