#ifndef MEMORYPOOL_H
#define MEMORYPOOL_H

#include <tuple>
#include <vector>

#include <Math/Factory.h>
#include <Math/Minimizer.h>
#include <TH1F.h>
#include <TVector3.h>

#include "Particle.h"

class alignas(64) MemoryPool {
protected:
  static thread_local Particle              *Detector_particleHitData_particle;
  static thread_local Particle              *Detector_reconstructUsingNonLinearMethod_particle;
  static thread_local std::vector<double>   *Detector_detect_detectedTimes;
  static thread_local std::vector<double>   *Detector_processReconstruction_hitTimes;
  static thread_local std::vector<double>   *Detector_processReconstruction_propagationLengths;
  static thread_local std::vector<double>   *Detector_reconstructUsingNonLinearMethod_reconstructedHitTimes;
  static thread_local ROOT::Math::Minimizer *Detector_reconstructUsingNonLinearMethod_minimizer;
  static thread_local std::vector<std::tuple<double, double, TVector3>> *Detector_particleHitData_hitData;

public:
  inline static Particle *getDetector_particleHitData_particle();

  inline static Particle *getDetector_reconstructUsingNonLinearMethod_particle();

  inline static std::vector<double> *getDetector_detect_detectedTimes(size_t size);

  inline static std::vector<double> *getDetector_processReconstruction_hitTimes(size_t size);

  inline static std::vector<double> *getDetector_processReconstruction_propagationLengths(size_t size);

  inline static std::vector<double> *getDetector_reconstructUsingNonLinearMethod_reconstructedHitTimes(size_t size);

  inline static ROOT::Math::Minimizer *getDetector_reconstructUsingNonLinearMethod_minimizer();

  inline static std::vector<std::tuple<double, double, TVector3>> *getDetector_particleHitData_hitData(size_t size);
};

thread_local Particle            *MemoryPool::Detector_particleHitData_particle                              = nullptr;
thread_local Particle            *MemoryPool::Detector_reconstructUsingNonLinearMethod_particle              = nullptr;
thread_local std::vector<double> *MemoryPool::Detector_detect_detectedTimes                                  = nullptr;
thread_local std::vector<double> *MemoryPool::Detector_processReconstruction_hitTimes                        = nullptr;
thread_local std::vector<double> *MemoryPool::Detector_processReconstruction_propagationLengths              = nullptr;
thread_local std::vector<double> *MemoryPool::Detector_reconstructUsingNonLinearMethod_reconstructedHitTimes = nullptr;
thread_local ROOT::Math::Minimizer *MemoryPool::Detector_reconstructUsingNonLinearMethod_minimizer           = nullptr;
thread_local std::vector<std::tuple<double, double, TVector3>> *MemoryPool::Detector_particleHitData_hitData = nullptr;

inline Particle *MemoryPool::getDetector_particleHitData_particle() {
  if (!Detector_particleHitData_particle) Detector_particleHitData_particle = new Particle();

  return Detector_particleHitData_particle;
};

inline Particle *MemoryPool::getDetector_reconstructUsingNonLinearMethod_particle() {
  if (!Detector_reconstructUsingNonLinearMethod_particle)
    Detector_reconstructUsingNonLinearMethod_particle = new Particle();

  return Detector_reconstructUsingNonLinearMethod_particle;
};

inline std::vector<double> *MemoryPool::getDetector_detect_detectedTimes(size_t size) {
  if (!Detector_detect_detectedTimes) Detector_detect_detectedTimes = new std::vector<double>(size);

  Detector_detect_detectedTimes->clear();
  return Detector_detect_detectedTimes;
};

inline std::vector<double> *MemoryPool::getDetector_processReconstruction_hitTimes(size_t size) {
  if (!Detector_processReconstruction_hitTimes) Detector_processReconstruction_hitTimes = new std::vector<double>(size);

  Detector_processReconstruction_hitTimes->clear();
  return Detector_processReconstruction_hitTimes;
};

inline std::vector<double> *MemoryPool::getDetector_processReconstruction_propagationLengths(size_t size) {
  if (!Detector_processReconstruction_propagationLengths)
    Detector_processReconstruction_propagationLengths = new std::vector<double>(size);

  Detector_processReconstruction_propagationLengths->clear();
  return Detector_processReconstruction_propagationLengths;
};

inline std::vector<double> *MemoryPool::getDetector_reconstructUsingNonLinearMethod_reconstructedHitTimes(size_t size) {
  if (!Detector_reconstructUsingNonLinearMethod_reconstructedHitTimes)
    Detector_reconstructUsingNonLinearMethod_reconstructedHitTimes = new std::vector<double>(size);

  Detector_reconstructUsingNonLinearMethod_reconstructedHitTimes->clear();
  return Detector_reconstructUsingNonLinearMethod_reconstructedHitTimes;
};

inline ROOT::Math::Minimizer *MemoryPool::getDetector_reconstructUsingNonLinearMethod_minimizer() {
  if (!Detector_reconstructUsingNonLinearMethod_minimizer)
    Detector_reconstructUsingNonLinearMethod_minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

  return Detector_reconstructUsingNonLinearMethod_minimizer;
};

inline std::vector<std::tuple<double, double, TVector3>> *MemoryPool::getDetector_particleHitData_hitData(size_t size) {
  if (!Detector_particleHitData_hitData)
    Detector_particleHitData_hitData = new std::vector<std::tuple<double, double, TVector3>>(size);

  Detector_particleHitData_hitData->clear();
  return Detector_particleHitData_hitData;
};

#endif // MEMORYPOOL_H