// Coordinate system definition
/********************************
 *          x -------·          *
 *                  /|          *
 *                 / |          *
 *                y  |          *
 *                   z          *
 *******************************/

#ifndef DETECTOR_H
#define DETECTOR_H

#include <string>
#include <vector>

#include <Math/Vector3D.h>
#include <TGraphErrors.h>
#include <TRandom3.h>

#include "Config.h"
#include "HitsData.h"
#include "Particle.h"
#include "ScintillatorCounters.h"
#include "ThreadPool.h"

class Detector {
public:
  /* BEGIN Constructor & Destructor */

  /**
   * @brief Constructor
   * @param scintillatorCounters Scintillator counters
   * @param B Magnetic field
   */
  Detector(const std::vector<ScintillatorCounters> &scintillatorCounters, const ROOT::Math::XYZVector &B);

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
  inline ROOT::Math::XYZVector getB() const;

  /* END Getters */

  /* BEGIN Methods */

  /**
   * @brief Calculate the hit times and positions of the particle in this detector
   * @param particle Incident particle
   * @param enableEnergyLoss Enable energy loss in the scintillator counters
   * @param enableEnergyLossFluctuation Enable energy loss fluctuation in the scintillator counters
   * (NOTE! This can be effective OLNY when enableEnergyLoss is true and Config::enableEnergyLossFluctuation is true)
   * @return Hits data of the particle in this detector
   */
  HitsData *particleHitData(
      const Particle &particleOrignal,
      const bool      enableEnergyLoss            = true,
      const bool      enableEnergyLossFluctuation = true,
      const HitsData *measuresData                = nullptr
  ) const;

  /**
   * @brief Plot the time difference between the hit times of the particle in this detector with and without energy loss
   * @param particle Incident particle
   * @param fileName Name of the file to save the plot
   */
  void plotDeltaTime(const Particle &particle, const std::string &fileName = "test.png") const;

  /**
   * @brief Plot the time difference between the hit times of the particle in this detector with and without energy loss
   * @param hitDataWithEnergyLoss Hits data of the particle in this detector with energy loss
   * @param hitDataWithoutEnergyLoss Hits data of the particle in this detector without energy loss
   * @param fileName Name of the file to save the plot
   */
  void plotDeltaTime(
      const HitsData    &hitDataWithEnergyLoss,
      const HitsData    &hitDataWithoutEnergyLoss,
      const std::string &fileName = "test.png"
  ) const;

  /**
   * @brief Detect the particle
   * @param particle Incident particle
   * @return Measures data of the particle in the scintillator counters
   */
  inline HitsData *measure(const Particle &particle) const;

  /**
   * @brief Detect the particle
   * @param hitsData Hits data of the particle in this detector
   * @return Measures data of the particle in the scintillator counters
   */
  HitsData *measure(const HitsData &hitsData) const;

  /**
   * @brief Reconstruct the particle's beta using the linear method
   * @param particle Incident particle
   * @return Reconstructed 1/beta of the particle
   */
  inline double reconstructUsingLinearMethod(const Particle &particle) const;

  /**
   * @brief Reconstruct the particle's beta using the linear method
   * @param measuresData Measures data of the particle in the scintillator counters
   * @return Reconstructed 1/beta of the particle
   */
  double reconstructUsingLinearMethod(const HitsData &measuresData) const;

  /**
   * @brief Plot the reconstructed data of the particle in this detector
   * @param particle Incident particle
   * @param fileName Name of the file to save the plot
   */
  void plotReconstructDataUsingLinearMethod(const Particle &particle, const std::string &fileName = "test.png") const;

  /**
   * @brief Calculate the distribution of the reconstructed 1/beta of the particle in this detector
   * @param particle Incident particle
   * @param nReconstructions Number of reconstructions
   * @param enableLinearMethod Enable reconstruction using the linear method
   * @param enablePlot Enable plot
   * @param fileName Name of the file to save the plot
   * @return Mean and standard deviation of the distribution of the reconstructed 1/beta of the particle
   */
  std::pair<double, double> distributionOfReconstruction(
      const Particle    &particle,
      const int          nReconstructions  = 10000,
      const bool         enableLinearMethod = false,
      const bool         enablePlot        = true,
      const std::string &fileName          = "test.png"
  ) const;

  /**
   * @brief Plot the difference between real and reconstructed 1/beta of the particle in this detectors
   * @param particle Incident particle
   * @param betaMin Minimum beta of the particle
   * @param betaMax Maximum beta of the particle
   * @param nPoints Number of points to plot (Infact, nPoints + 1 points will be plotted)
   * @param enableLinearMethod Enable reconstruction using the linear method
   * @param fileName Name of the file to save the plot
   * @param enablePlot Enable plot
   * @return Graph of the difference between real and reconstructed 1/beta of the particle
   */
  TGraphErrors *deltaBetaReciprocal(
      Particle           particle,
      const double       betaMin            = 0.4,
      const double       betaMax            = 0.9,
      const int          nPoints            = 5,
      const bool         enableLinearMethod = false,
      const bool         enablePlot         = false,
      const std::string &fileName           = "test.png"
  ) const;

  /**
   * @brief Reconstruct the particle's beta using the non-linear method
   * @param particle Incident particle
   * @return Reconstructed 1/beta of the particle
   */
  inline double reconstructUsingNonLinearMethod(const Particle &particle) const;

  /**
   * @brief Reconstruct the particle's beta using the non-linear method
   * @param particle Incident particle
   * @param measuresData Measures data of the particle in the scintillator counters
   * @return Reconstructed 1/beta of the particle
   */
  double reconstructUsingNonLinearMethod(const Particle &particleOrignal, const HitsData &measuresDataOrignal) const;

  /**
   * @brief Plot the trajectory of the particle in this detector
   * @param particle Incident particle
   * @param fileName Name of the file to save the plot
   */
  void plotParticleTrajectory(Particle particle, const std::string &fileName = "test.png") const;

  /* END Methods */

protected:
  std::vector<ScintillatorCounters> dScintillatorCounters; // Scintillator counters
  ROOT::Math::XYZVector             dB;                    // Magnetic field
};

inline double Detector::getMinZ() const { return this->dScintillatorCounters.front().getLocation(); }

inline std::vector<ScintillatorCounters> Detector::getScintillatorCounters() const {
  return this->dScintillatorCounters;
}

inline ROOT::Math::XYZVector Detector::getB() const { return this->dB; }

inline HitsData *Detector::measure(const Particle &particle) const {
  return this->measure(*this->particleHitData(particle, true));
}

inline double Detector::reconstructUsingLinearMethod(const Particle &particle) const {
  const HitsData *hitsData = this->particleHitData(particle, true);
  return this->reconstructUsingLinearMethod(*this->measure(*hitsData));
}

inline double Detector::reconstructUsingNonLinearMethod(const Particle &particle) const {
  const HitsData *hitsData = this->particleHitData(particle, true);
  return this->reconstructUsingNonLinearMethod(particle, *this->measure(*hitsData));
}

#endif /* DETECTOR_H */
