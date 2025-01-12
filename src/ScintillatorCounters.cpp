#include "ScintillatorCounters.h"

ScintillatorCounters::ScintillatorCounters(double location, bool direction)
    : location(location), direction(direction) {}

ScintillatorCounters::~ScintillatorCounters() {}

double ScintillatorCounters::getLocation() const
{
    return location;
}

bool ScintillatorCounters::getDirection() const
{
    return direction;
}