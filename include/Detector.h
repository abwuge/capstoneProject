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

#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <TGraphErrors.h>
#include <TRandom3.h>
#include <TVector3.h>

#include "Config.h"
#include "Particle.h"
#include "ScintillatorCounters.h"
#include "ThreadPool.h"

class Detector {
protected:
  std::vector<ScintillatorCounters> scintillatorCounters; // Scintillator counters
  TVector3                          B;                    // Magnetic field

  /* BEGIN Methods */

  /**
   * @brief Process the reconstruction of the particle
   * @param particle Incident particle
   * @param enableLinerMethod Enable the linear method
   * @param betaReciprocalReal Real 1/beta of the particle
   * @param results Results of the reconstruction
   * @param index Index of the result
   */
  void processReconstruction(
      const Particle      &particle,
      const bool           enableLinerMethod,
      const double         betaReciprocalReal,
      std::vector<double> &results,
      int                  index
  ) const;

  /* END Methods */

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
   * @brief Calculate the hit times and positions of the particle in this detector
   * @param particle Incident particle
   * @param enableEnergyLoss Enable energy loss in the scintillator counters
   * @param enableEnergyLossFluctuation Enable energy loss fluctuation in the scintillator counters (NOTE! This can be
   * effective OLNY when enableEnergyLoss is true and Config::enableEnergyLossFluctuation is true)
   * @return Hit times (in ns), propagation lengths (in cm), and hit position (in cm) of the particle in this detector
   */
  std::vector<std::tuple<double, double, TVector3>> *particleHitData(
      const Particle &particleOrignal,
      const bool      enableEnergyLoss            = true,
      const bool      enableEnergyLossFluctuation = true
  ) const;

  /**
   * @brief Plot the time difference between the hit times of the particle in this detector with and without energy loss
   *
   * This method uses TCanvas! So you need be CAREFUL where you call this method or remember to use cd() method in
   * TCanvas! Otherwise, you may get a blank plot with your TCanvas!
   * @param particle Incident particle
   * @param fileName Name of the file to save the plot
   */
  void plotDeltaTime(const Particle &particle, const std::string &fileName = "test.png") const;

  /**
   * @brief Plot the time difference between the hit times of the particle in this detector with and without energy loss
   *
   * This method uses TCanvas! So you need be CAREFUL where you call this method or remember to use cd() method in
   * TCanvas! Otherwise, you may get a blank plot with your TCanvas!
   * @param hitDataWithEnergyLoss Hit times (in ns), propagation lengths (in cm), and hit position (in cm) of the
   * particle in this detector with energy loss
   * @param hitDataWithoutEnergyLoss Hit times (in ns), propagation lengths (in cm), and hit position (in cm) of the
   * particle in this detector without energy loss
   * @param fileName Name of the file to save the plot
   */
  void plotDeltaTime(
      const std::vector<std::tuple<double, double, TVector3>> &hitDataWithEnergyLoss,
      const std::vector<std::tuple<double, double, TVector3>> &hitDataWithoutEnergyLoss,
      const std::string                                       &fileName = "test.png"
  ) const;

  /**
   * @brief Detect the particle
   *
   * For better performance, use hitTimes instead of particle
   * @param particle Incident particle
   * @return Detected times of the particle in the scintillator counters
   */
  std::vector<double> *detect(const Particle &particle) const;

  /**
   * @brief Detect the particle
   * @param hitTimes Hit times of the particle
   * @return Detected times of the particle in the scintillator counters
   */
  std::vector<double> *detect(const std::vector<double> &hitTimes) const;

  /**
   * @brief Reconstruct the particle's beta using the linear method
   * @param particle Incident particle
   * @return Reconstructed 1/beta of the particle
   */
  double reconstructUsingLinearMethod(const Particle &particle) const;

  /**
   * @brief Reconstruct the particle's beta using the linear method
   * @param detectedTimes Detected times of the particle in the scintillator counters
   * @param propagationLengths Propagation lengths of the particle in this decetor
   * @return Reconstructed 1/beta of the particle
   */
  double reconstructUsingLinearMethod(
      const std::vector<double> &detectedTimes,
      const std::vector<double> &propagationLengths
  ) const;

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
   * @param enableLinerMethod Enable reconstruction using the linear method
   * @param enablePlot Enable plot
   * @param fileName Name of the file to save the plot
   * @return Mean and standard deviation of the distribution of the reconstructed 1/beta of the particle
   */
  std::pair<double, double> distributionOfReconstruction(
      const Particle    &particle,
      const int          nReconstructions  = 10000,
      const bool         enableLinerMethod = false,
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
  double reconstructUsingNonLinearMethod(const Particle &particle) const;

  /**
   * @brief Reconstruct the particle's beta using the non-linear method
   * @param particle Incident particle
   * @param detectedTimes Detected times of the particle in the scintillator counters
   * @param propagationLengths Propagation lengths of the particle in this decetor
   * @return Reconstructed 1/beta of the particle
   */
  double reconstructUsingNonLinearMethod(
      const Particle            &particleOrignal,
      const std::vector<double> &detectedTimes,
      const std::vector<double> &propagationLengths
  ) const;

  /* END Methods */
};

inline double Detector::getMinZ() const { return this->scintillatorCounters.front().getLocation(); }

inline std::vector<ScintillatorCounters> Detector::getScintillatorCounters() const {
  return this->scintillatorCounters;
}

inline TVector3 Detector::getB() const { return this->B; }

#endif /* DETECTOR_H */
