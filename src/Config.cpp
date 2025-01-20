#include "Config.h"

// Infact, we can use TRandom3 *Detector::Random = new TRandom3(configEnableFixedSeed); to implement the following code
// But the following code is more readable
#if configEnableFixedSeed
TRandom3 *Config::Random = new TRandom3(1); // Fixed seed 1
#else
TRandom3 *Config::Random = new TRandom3(0); // 0 means time-based seed
#endif