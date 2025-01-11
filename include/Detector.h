// Coordinate system definition
/*******************************
*           ·------- x         *
*          /|                  *
*         / |                  *
*        y  |                  *
*           z                  *
*******************************/

#ifndef DETECTOR_H
#define DETECTOR_H

#include <vector>
#include <TVector3.h>

class Detector
{
private:
    std::vector<double> scintillatorCountersLocation; // z-coordinates of the scintillator counters in cm (assuming they are infinite planes in the x and y directions)
    TVector3 B;                                       // Magnetic field

public:
    // Constructor & Destructor
    Detector(const std::vector<double> &scintillatorCountersLocation, const TVector3 &B);
    ~Detector();

    // Getters
    double getMinZ() const;                                      // Minimum z-coordinate of the scintillator counters in cm
    std::vector<double> getScintillatorCountersLocation() const; // z-coordinates of the scintillator counters in cm (assuming they are infinite planes in the x and y directions)
    TVector3 getB() const;                                       // Magnetic field
};

#endif /* DETECTOR_H */