#ifndef SCINTILLATORCOUNTERS_H
#define SCINTILLATORCOUNTERS_H

#include "Material.h"
#include "Particle.h"

class ScintillatorCounters
{
private:
    double location;   // z-coordinate of the scintillator counter in cm (assuming it is an infinite plane in the x and y directions)
    bool direction;    // Indicates whether the detector detects x (true) or y (false) coordinates
    double thickness;  // Thickness of the scintillator counter in cm
    Material material; // Material of the scintillator counter

public:
    /* BEGIN Constructor & Destructor */

    /**
     * @brief Constructor
     * @param location z-coordinate of the scintillator counter in cm (assuming it is an infinite plane in the x and y directions)
     * @param direction Indicates whether the detector detects x (true) or y (false) coordinates
     * @param thickness Thickness of the scintillator counter in cm
     * @param material Material of the scintillator counter
     */
    ScintillatorCounters(const double location, const bool direction, const double thickness, const Material &material);
    
    ~ScintillatorCounters();

    /* END Constructor & Destructor */

    /* BEGIN Getters */

    /**
     * @brief Get the z-coordinate of the scintillator counter
     */
    double getLocation() const;

    /**
     * @brief Get whether the detector detects x (true) or y (false) coordinates
     */
    bool getDirection() const;

    /**
     * @brief Get the thickness of the scintillator counter
     */
    double getThickness() const;

    /**
     * @brief Get the material of the scintillator counter
     */
    Material getMaterial() const;

    /* END Getters */

    /* BEGIN Methods */

    /**
     * @brief Calculate the energy loss of the particle in the scintillator counter
     * @param particle Particle
     * @return Energy loss of the particle in the scintillator counter in MeV
     */
    double energyLoss(const Particle &particle) const;

    /* END Methods */
};

#endif /* SCINTILLATORCOUNTERS_H */
