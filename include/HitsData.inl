#include "HitsData.h"

inline HitsData::HitsData(size_t size) {
  hitsData.reserve(size);
  hitsTime.reserve(size);
  hitsLength.reserve(size);
  hitsEnergyLoss.reserve(size);
  hitsPosition.reserve(size);
}

inline const HitData &HitsData::at(size_t index) const { return hitsData.at(index); }

inline const std::vector<double> &HitsData::getHitsTime() const { return hitsTime; };

inline const std::vector<double> &HitsData::getHitsLength() const { return hitsLength; };

inline const std::vector<double> &HitsData::getHitsEnergyLoss() const { return hitsEnergyLoss; };

inline const std::vector<ROOT::Math::XYZVector> &HitsData::getHitsPosition() const { return hitsPosition; };

inline size_t HitsData::size() const { return hitsData.size(); }

inline void HitsData::push_back(const HitData &hitData) {
  hitsData.push_back(hitData);
  hitsTime.push_back(hitData.hitTime);
  hitsLength.push_back(hitData.hitLength);
  hitsEnergyLoss.push_back(hitData.hitEnergyLoss);
  hitsPosition.push_back(hitData.hitPosition);
}

inline void HitsData::set(size_t index, const HitData &hitData) {
  hitsData.at(index)       = hitData;
  hitsTime.at(index)       = hitData.hitTime;
  hitsLength.at(index)     = hitData.hitLength;
  hitsEnergyLoss.at(index) = hitData.hitEnergyLoss;
  hitsPosition.at(index)   = hitData.hitPosition;
}

inline void HitsData::setHitTime(size_t index, double hitTime) {
  hitsTime.at(index)         = hitTime;
  hitsData.at(index).hitTime = hitTime;
}

inline void HitsData::setHitLength(size_t index, double hitLength) {
  hitsLength.at(index)         = hitLength;
  hitsData.at(index).hitLength = hitLength;
}

inline void HitsData::setHitEnergyLoss(size_t index, double hitEnergyLoss) {
  hitsEnergyLoss.at(index)         = hitEnergyLoss;
  hitsData.at(index).hitEnergyLoss = hitEnergyLoss;
}

inline void HitsData::setHitPosition(size_t index, const ROOT::Math::XYZVector &hitPosition) {
  hitsPosition.at(index)         = hitPosition;
  hitsData.at(index).hitPosition = hitPosition;
}

inline void HitsData::push_back(
    const double                 hitTime,
    const double                 hitLength,
    const double                 hitEnergyLoss,
    const ROOT::Math::XYZVector &hitPosition
) {
  hitsData.push_back({hitTime, hitLength, hitEnergyLoss, hitPosition});
  hitsTime.push_back(hitTime);
  hitsLength.push_back(hitLength);
  hitsEnergyLoss.push_back(hitEnergyLoss);
  hitsPosition.push_back(hitPosition);
}

inline void HitsData::clear() {
  hitsData.clear();
  hitsTime.clear();
  hitsLength.clear();
  hitsEnergyLoss.clear();
  hitsPosition.clear();
}