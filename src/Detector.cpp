#include "Detector.h"
#include <algorithm>

Detector::Detector(const std::vector<ScintillatorCounters> &scintillatorCounters, const TVector3 &B)
{
    this->scintillatorCounters = scintillatorCounters;
    this->B = B;
    std::sort(this->scintillatorCounters.begin(), this->scintillatorCounters.end(),
              [](const ScintillatorCounters &a, const ScintillatorCounters &b)
              { return a.getLocation() < b.getLocation(); });
}

Detector::~Detector() {}

double Detector::getMinZ() const
{
    return this->scintillatorCounters.front().getLocation();
}

std::vector<ScintillatorCounters> Detector::getScintillatorCounters() const
{
    return this->scintillatorCounters;
}

TVector3 Detector::getB() const
{
    return this->B;
}
