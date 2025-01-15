#ifndef SCINTILLATORCOUNTERS_H
#define SCINTILLATORCOUNTERS_H

#include <string>

#include "config.h"
#include "Material.h"
#include "Particle.h"

class ScintillatorCounters
{
private:
    double location;       // z-coordinate of the scintillator counter in cm (assuming it is an infinite plane in the x and y directions)
    bool direction;        // Indicates whether the detector detects x (true) or y (false) coordinates
    double thickness;      // Thickness of the scintillator counter in cm
    double timeResolution; // Time resolution of the scintillator counter in ns
    Material material;     // Material of the scintillator counter

public:
    /* BEGIN Constructor & Destructor */

    /**
     * @brief Constructor
     * @param location z-coordinate of the scintillator counter in cm (assuming it is an infinite plane in the x and y directions)
     * @param direction Indicates whether the detector detects x (true) or y (false) coordinates
     * @param thickness Thickness of the scintillator counter in cm
     * @param material Material of the scintillator counter
     */
    ScintillatorCounters(const double location, const bool direction, const double thickness, const double timeResolution, const Material &material);

    ~ScintillatorCounters();

    /* END Constructor & Destructor */

    /* BEGIN Getters */

    /**
     * @brief Get the z-coordinate of the scintillator counter
     */
    inline double getLocation() const;

    /**
     * @brief Get whether the detector detects x (true) or y (false) coordinates
     */
    inline bool getDirection() const;

    /**
     * @brief Get the thickness of the scintillator counter
     */
    inline double getThickness() const;

    /**
     * @brief Get the time resolution of the scintillator counter
     */
    inline double getTimeResolution() const;

    /**
     * @brief Get the material of the scintillator counter
     */
    inline Material getMaterial() const;

    /* END Getters */

    /* BEGIN Methods */

    /**
     * @brief Calculate the energy loss of the particle in the scintillator counter
     * @param particle Particle
     * @return Energy loss of the particle in the scintillator counter in MeV
     */
    double energyLoss(const Particle &particle) const;

    /**
     * @brief Plot the energy loss of the particle in the scintillator counter using the Bethe-Bloch formula
     *
     * This method uses TCanvas! So you need be CAREFUL where you call this method or remember to use cd() method in TCanvas! Otherwise, you may get a blank plot with your TCanvas!
     * @param particle Particle
     * @param betaGammaMin Minimum beta * gamma of the particle
     * @param betaGammaMax Maximum beta * gamma of the particle
     * @param nPoints Number of points to plot (Infact, nPoints + 1 points will be plotted)
     * @param enableKineticEnergy Enable kinetic energy in the plot
     * @param fileName Name of the file to save the plot
     */
    void plotEnergyLoss(Particle particle, const double betaGammaMin = 0.1, const double betaGammaMax = 1000, const int nPoints = 1000, const bool enableKineticEnergy = true, const std::string &fileName = "test.png") const;

    /* END Methods */
};

inline double ScintillatorCounters::getLocation() const { return this->location; }

inline bool ScintillatorCounters::getDirection() const { return this->direction; }

inline double ScintillatorCounters::getThickness() const { return this->thickness; }

inline double ScintillatorCounters::getTimeResolution() const { return this->timeResolution; }

inline Material ScintillatorCounters::getMaterial() const { return this->material; }

#endif /* SCINTILLATORCOUNTERS_H */
