#include "Detector.h"
#include <algorithm>

Detector::Detector(const std::vector<double> &scintillatorCountersLocation, const TVector3 &B)
{
    this->scintillatorCountersLocation = scintillatorCountersLocation;
    this->B = B;
    std::sort(this->scintillatorCountersLocation.begin(), this->scintillatorCountersLocation.end());
}

Detector::~Detector() {}

double Detector::getMinZ() const
{
    return this->scintillatorCountersLocation.front();
}

std::vector<double> Detector::getScintillatorCountersLocation() const
{
    return this->scintillatorCountersLocation;
}

TVector3 Detector::getB() const
{
    return this->B;
}