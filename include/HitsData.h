#ifndef HITDATA_H
#define HITDATA_H

#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

#include <TVector3.h>

struct HitData {
  double   hitTime;
  double   hitLength;
  double   hitEnergyLoss;
  TVector3 hitPosition;
};

class HitsData {
public:
  /* BEGIN Constructor & Destructor */

  /**
   * @brief Constructor
   * @param size Size of the hits data
   */
  inline HitsData(size_t size = 4);

  ~HitsData() {}

  /* END Constructor & Destructor */

  /* BEGIN Getters */

  /**
   * @brief Get the hit data at the index
   * @param index Index
   * @return Hit data
   */
  inline const HitData &at(size_t index) const;

  /**
   * @brief Get hits time
   */
  inline const std::vector<double> &getHitsTime() const;

  /**
   * @brief Get propagation lengths
   */
  inline const std::vector<double> &getHitsLength() const;

  /**
   * @brief Get energy loss
   */
  inline const std::vector<double> &getHitsEnergyLoss() const;

  /**
   * @brief Get hit positions
   */
  inline const std::vector<TVector3> &getHitsPosition() const;

  /**
   * @brief Get the size of the hits data
   */
  inline size_t size() const;

  /* END Getters */

  /* BEGIN Setters*/

  /**
   * @brief Set the hit data at the index
   * @param index Index
   * @param hitData Hit data
   */
  inline void set(size_t index, const HitData &hitData);

  /**
   * @brief Set the hit time at the index
   * @param index Index
   * @param hitTime Hit time (in ns)
   */
  inline void setHitTime(size_t index, double hitTime);

  /**
   * @brief Set the propagation length at the index
   * @param index Index
   * @param hitLength Propagation length (in cm)
   */
  inline void setHitLength(size_t index, double hitLength);

  /**
   * @brief Set the energy loss at the index
   * @param index Index
   * @param hitEnergyLoss Energy loss (in MeV)
   */
  inline void setHitEnergyLoss(size_t index, double hitEnergyLoss);

  /**
   * @brief Set the hit position at the index
   * @param index Index
   * @param hitPosition Hit position (in cm)
   */
  inline void setHitPosition(size_t index, const TVector3 &hitPosition);

  /* END Setters */

  /* BEGIN Methods */

  /**
   * @brief Add a hit data to the hits data
   * @param hitData Hit data
   */
  inline void push_back(const HitData &hitData);

  /**
   * @brief Add a hit data to the hits data
   * @param hitTime Hit time (in ns)
   * @param hitLength Propagation length (in cm)
   * @param hitEnergyLoss Energy loss (in MeV)
   * @param hitPosition Hit position (in cm)
   */
  inline void
  push_back(const double hitTime, const double hitLength, const double hitEnergyLoss, const TVector3 &hitPosition);

  /**
   * @brief Clear the hits data
   */
  inline void clear();

  /* END Methods */

protected:
  std::vector<HitData>  hitsData;       // Hit data
  std::vector<double>   hitsTime;       // Hit times (in ns)
  std::vector<double>   hitsLength;     // Propagation lengths (in cm)
  std::vector<double>   hitsEnergyLoss; // Energy loss (in MeV)
  std::vector<TVector3> hitsPosition;   // Hit position (in cm)
};

#include "HitsData.inl"

#endif /* HITDATA_H */