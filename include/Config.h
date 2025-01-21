#ifndef CONFIG_H
#define CONFIG_H

// random related
#define configEnableFixedSeed 1

// notice related
#define configEnableWarning 0
#define configEnableWarningAll 0

// test related
#define configEnableTest 0
#define configEnableDebug 0
#define configEnableDebugAtCoding 1

// todo related
#define configTodo 1

// influence related
#define configEnableTimeResolution 1
#define configUseBetheBloch 0
#define configEnableEnergyLossFluctuation 1

// multi-threading related
#define configEnableMultiThreading 1
#define configEnableMultiThreadingAnywhere 0

#include <TRandom3.h>
class Config
{
public:
    static TRandom3 *Random;
};

#endif /* CONFIG_H */