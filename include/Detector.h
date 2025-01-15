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
#include <string>

#include <TVector3.h>
#include <TRandom3.h>

#include "config.h"
#include "Particle.h"
#include "ScintillatorCounters.h"

class Detector
{
private:
    static TRandom3 *Random;                                 // Random number generator
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

    /**
     * @brief Plot the time difference between the hit times of the particle in the scintillator counters with and without energy loss
     *
     * This method uses TCanvas! So you need be CAREFUL where you call this method or remember to use cd() method in TCanvas! Otherwise, you may get a blank plot with your TCanvas!
     * @param particle Incident particle
     * @param fileName Name of the file to save the plot
     */
    void plotDeltaTime(const Particle &particle, const std::string &fileName = "test.png") const;

    /**
     * @brief Plot the time difference between the hit times of the particle in the scintillator counters with and without energy loss
     *
     * This method uses TCanvas! So you need be CAREFUL where you call this method or remember to use cd() method in TCanvas! Otherwise, you may get a blank plot with your TCanvas!
     * @param hitDataWithEnergyLoss Hit time (in ns), propagation length (in cm), and hit position (in cm) of the particle in the scintillator counters with energy loss
     * @param hitDataWithoutEnergyLoss Hit time (in ns), propagation length (in cm), and hit position (in cm) of the particle in the scintillator counters without energy loss
     * @param fileName Name of the file to save the plot
     */
    void plotDeltaTime(const std::vector<std::tuple<double, double, TVector3>> &hitDataWithEnergyLoss, const std::vector<std::tuple<double, double, TVector3>> &hitDataWithoutEnergyLoss, const std::string &fileName = "test.png") const;

    /**
     * @brief Detect the particle in the scintillator counters
     * 
     * For better performance, use hitTimes instead of particle
     * @param particle Incident particle
     * @return Hit times of the particle in the scintillator counters
     */
    std::vector<double> detect(const Particle &particle) const;

    /**
     * @brief Detect the particle in the scintillator counters
     * @param hitTimes Hit times of the particle
     * @return Detected times of the particle in the scintillator counters
     */
    std::vector<double> detect(const std::vector<double> &hitTimes) const;

    /**
     * @brief Reconstruct the particle's beta using the linear method
     * @param particle Incident particle
     * @return Reconstructed beta of the particle
     */
    double reconstructUsingLinearMethod(const Particle &particle) const;

    /**
     * @brief Reconstruct the particle's beta using the linear method
     * @param hitTimes Hit times of the particle
     * @return Reconstructed beta of the particle
     */
    double reconstructUsingLinearMethod(const std::vector<double> &hitTimes) const;

    /* END Methods */
};

inline double Detector::getMinZ() const { return this->scintillatorCounters.front().getLocation(); }

inline std::vector<ScintillatorCounters> Detector::getScintillatorCounters() const { return this->scintillatorCounters; }

inline TVector3 Detector::getB() const { return this->B; }

#endif /* DETECTOR_H */
