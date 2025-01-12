#include "Detector.h"
#include <algorithm>

Detector::Detector(const std::vector<ScintillatorCounters> &scintillatorCounters, const TVector3 &B)
    : scintillatorCounters(scintillatorCounters), B(B)
{
    std::sort(this->scintillatorCounters.begin(), this->scintillatorCounters.end(),
              [](const ScintillatorCounters &a, const ScintillatorCounters &b)
              { return a.getLocation() < b.getLocation(); });
}

Detector::~Detector() {}

double Detector::getMinZ() const
{
    return scintillatorCounters.front().getLocation();
}

std::vector<ScintillatorCounters> Detector::getScintillatorCounters() const
{
    return scintillatorCounters;
}

TVector3 Detector::getB() const
{
    return B;
}
