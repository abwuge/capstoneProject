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

#include "ScintillatorCounters.h"

class Detector
{
private:
    std::vector<ScintillatorCounters> scintillatorCounters; // Scintillator counters
    TVector3 B;                                             // Magnetic field

public:
    /* BEGIN Constructor & Destructor */
    Detector(const std::vector<ScintillatorCounters> &scintillatorCounters, const TVector3 &B);
    ~Detector();
    /* END Constructor & Destructor */

    /* BEGIN Getters */
    // returns the minimum z-coordinate of the scintillator counters in cm
    double getMinZ() const;

    // returns the scintillator counters
    std::vector<ScintillatorCounters> getScintillatorCounters() const;

    // returns the magnetic field
    TVector3 getB() const;
    /* END Getters */
};

#endif /* DETECTOR_H */
