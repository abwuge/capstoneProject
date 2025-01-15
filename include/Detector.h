// Coordinate system definition
/*******************************
 *           ·------- x         *
 *          /|                  *
 *         / |                  *
 *        y  |                  *
 *           z                  *
 *******************************/

#ifndef DETECTOR_H
#define DETECTOR_H

#include <tuple>
#include <vector>

#include <TVector3.h>

#include "config.h"
#include "Particle.h"
#include "ScintillatorCounters.h"

class Detector
{
private:
    std::vector<ScintillatorCounters> scintillatorCounters; // Scintillator counters
    TVector3 B;                                             // Magnetic field

public:
    /* BEGIN Constructor & Destructor */

    /**
     * @brief Constructor
     * @param scintillatorCounters Scintillator counters
     * @param B Magnetic field
     */
    Detector(const std::vector<ScintillatorCounters> &scintillatorCounters, const TVector3 &B);

    ~Detector();

    /* END Constructor & Destructor */

    /* BEGIN Getters */

    /**
     * @brief Get the minimum z-coordinate of the scintillator counters
     */
    inline double getMinZ() const;

    /**
     * @brief Get the scintillator counters
     */
    inline std::vector<ScintillatorCounters> getScintillatorCounters() const;

    /**
     * @brief Get the magnetic field
     */
    inline TVector3 getB() const;

    /* END Getters */

    /* BEGIN Methods */

    /**
     * @brief Calculate the cyclotron radius of the particle in the magnetic field
     * @param particle Particle
     * @return Cyclotron radius of the particle in the magnetic field
     */
    double particleCyclotronRadius(const Particle &particle) const;

    /**
     * @brief Calculate the direction of the particle in the magnetic field
     * @param particle Particle
     * @return Direction of the particle in the magnetic field
     */
    TVector3 particleCyclotronDirection(const Particle &particle) const;

    /**
     * @brief Calculate the hit time and positions of the particle in the scintillator counters
     * @param particle Incident particle
     * @param enableEnergyLoss Enable energy loss in the scintillator counters
     * @return Hit time (in ns), propagation length (in cm), and hit position (in cm) of the particle in the scintillator counters
     */
    std::vector<std::tuple<double, double, TVector3>> particleHitData(Particle particle, const bool enableEnergyLoss = true) const;

    /* END Methods */
};

inline double Detector::getMinZ() const { return this->scintillatorCounters.front().getLocation(); }

inline std::vector<ScintillatorCounters> Detector::getScintillatorCounters() const { return this->scintillatorCounters; }

inline TVector3 Detector::getB() const { return this->B; }

#endif /* DETECTOR_H */
