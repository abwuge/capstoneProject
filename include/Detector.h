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

#include <vector>
#include <TVector3.h>

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
    double getMinZ() const;

    /**
     * @brief Get the scintillator counters
     */
    std::vector<ScintillatorCounters> getScintillatorCounters() const;

    /**
     * @brief Get the magnetic field
     */
    TVector3 getB() const;

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
     * @brief Calculate the hit positions of the particle in the scintillator counters
     * @param particle Particle
     * @return Hit positions of the particle in the scintillator counters
     */
    std::vector<TVector3> particleHitPositions(const Particle &particle) const;

    /* END Methods */
};

#endif /* DETECTOR_H */
