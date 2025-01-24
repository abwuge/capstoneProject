#ifndef MEMORYPOOL_H
#define MEMORYPOOL_H

#include <Math/Factory.h>
#include <Math/Minimizer.h>

#include "HitsData.h"
#include "Particle.h"

class alignas(64) MemoryPool {
protected:
  static thread_local Particle              *Detector_particleHitData_particle;
  static thread_local Particle              *Detector_reconstructUsingNonLinearMethod_particle;
  static thread_local HitsData              *Detector_particleHitData_hitsData;
  static thread_local HitsData              *Detector_detect_detectedsData;
  static thread_local HitsData              *Detector_distributionOfReconstruction_hitsData;
  static thread_local HitsData              *Detector_reconstructUsingNonLinearMethod_measuresData;
  static thread_local ROOT::Math::Minimizer *Detector_reconstructUsingNonLinearMethod_minimizer;

public:
  inline static Particle *getDetector_particleHitData_particle();

  inline static Particle *getDetector_reconstructUsingNonLinearMethod_particle();

  inline static HitsData *getDetector_particleHitData_hitsData(size_t size);

  inline static HitsData *getDetector_detect_detectedsData(size_t size);

  inline static HitsData *getDetector_distributionOfReconstruction_hitsData(size_t size);

  inline static HitsData *getDetector_reconstructUsingNonLinearMethod_measuresData(size_t size);

  inline static ROOT::Math::Minimizer *getDetector_reconstructUsingNonLinearMethod_minimizer();
};

thread_local Particle              *MemoryPool::Detector_particleHitData_particle                     = nullptr;
thread_local Particle              *MemoryPool::Detector_reconstructUsingNonLinearMethod_particle     = nullptr;
thread_local HitsData              *MemoryPool::Detector_particleHitData_hitsData                     = nullptr;
thread_local HitsData              *MemoryPool::Detector_detect_detectedsData                         = nullptr;
thread_local HitsData              *MemoryPool::Detector_distributionOfReconstruction_hitsData        = nullptr;
thread_local HitsData              *MemoryPool::Detector_reconstructUsingNonLinearMethod_measuresData = nullptr;
thread_local ROOT::Math::Minimizer *MemoryPool::Detector_reconstructUsingNonLinearMethod_minimizer    = nullptr;

inline Particle *MemoryPool::getDetector_particleHitData_particle() {
  if (!Detector_particleHitData_particle) Detector_particleHitData_particle = new Particle();

  return Detector_particleHitData_particle;
};

inline Particle *MemoryPool::getDetector_reconstructUsingNonLinearMethod_particle() {
  if (!Detector_reconstructUsingNonLinearMethod_particle)
    Detector_reconstructUsingNonLinearMethod_particle = new Particle();

  return Detector_reconstructUsingNonLinearMethod_particle;
};

inline HitsData *MemoryPool::getDetector_detect_detectedsData(size_t size) {
  if (!Detector_detect_detectedsData) Detector_detect_detectedsData = new HitsData(size);

  Detector_detect_detectedsData->clear();
  return Detector_detect_detectedsData;
};

inline HitsData *MemoryPool::getDetector_distributionOfReconstruction_hitsData(size_t size) {
  if (!Detector_distributionOfReconstruction_hitsData)
    Detector_distributionOfReconstruction_hitsData = new HitsData(size);

  Detector_distributionOfReconstruction_hitsData->clear();
  return Detector_distributionOfReconstruction_hitsData;
};

inline HitsData *MemoryPool::getDetector_reconstructUsingNonLinearMethod_measuresData(size_t size) {
  if (!Detector_reconstructUsingNonLinearMethod_measuresData)
    Detector_reconstructUsingNonLinearMethod_measuresData = new HitsData(size);

  Detector_reconstructUsingNonLinearMethod_measuresData->clear();
  return Detector_reconstructUsingNonLinearMethod_measuresData;
};

inline ROOT::Math::Minimizer *MemoryPool::getDetector_reconstructUsingNonLinearMethod_minimizer() {
  if (!Detector_reconstructUsingNonLinearMethod_minimizer)
    Detector_reconstructUsingNonLinearMethod_minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

  return Detector_reconstructUsingNonLinearMethod_minimizer;
};

inline HitsData *MemoryPool::getDetector_particleHitData_hitsData(size_t size) {
  if (!Detector_particleHitData_hitsData) Detector_particleHitData_hitsData = new HitsData(size);

  Detector_particleHitData_hitsData->clear();
  return Detector_particleHitData_hitsData;
};

#endif // MEMORYPOOL_H