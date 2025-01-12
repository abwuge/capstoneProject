#include "ScintillatorCounters.h"

ScintillatorCounters::ScintillatorCounters(double location, bool direction)
{
    this->location = location;
    this->direction = direction;
}

ScintillatorCounters::~ScintillatorCounters() {}

double ScintillatorCounters::getLocation() const
{
    return this->location;
}

bool ScintillatorCounters::getDirection() const
{
    return this->direction;
}