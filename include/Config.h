#ifndef CONFIG_H
#define CONFIG_H

#include <TRandom3.h>
class Config {
public:
  static inline Config   *config = nullptr;               // Instance of the Config class
  static inline TRandom3 *random = nullptr;               // random number generator

  static inline bool enableFixedSeed = true;              // Enable fixed seed for random number generator

  static inline bool enableTest  = false;                 // Enable test
  static inline bool enableDebug = false;                 // Enable debug

  static inline bool enableWarning        = false;        // Enable warning
  static inline bool enableWarningAll     = false;        // Enable warning for all
  static inline bool enableWarningAtDebug = true;         // Enable warning at debug

  static inline bool todo = true;                         // Enable todo

  static inline bool enableTimeResolution        = true;  // Enable time resolution
  static inline bool useBetheBloch               = false; // Use Bethe-Bloch formula
  static inline bool enableEnergyLossFluctuation = true;  // Enable energy loss fluctuation
  static inline bool useLandau                   = true;  // Use Landau distribution
  static inline bool enableEnergyLossDetection   = true;  // Enable energy loss detection
  static inline bool useEnergyLoss2Rec           = true;  // use energy loss to reconstruct 1/ beta
  static inline bool useRealEnergyLoss2Rec       = true; // Use real energy loss

  static inline bool enableMultiThreading = true;         // Enable multi-threading

  /**
   * @brief Get the instance of the Config class
   */
  static inline Config *getInstance();

  /**
   * @brief Get the random number generator
   */
  static inline TRandom3 *getRandom();

protected:
  inline Config();

  Config(const Config &) = delete;

  Config &operator=(const Config &) = delete;
};

#include "Config.inl"

#endif /* CONFIG_H */