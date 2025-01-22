#include "Config.h"

inline Config *Config::getInstance() {
  static Config instance;
  return &instance;
}

inline Config::Config() {}

inline TRandom3 *Config::getRandom() {
  if (!random) random = new TRandom3(enableFixedSeed ? 1 : 0);

  return random;
}