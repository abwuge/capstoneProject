#ifndef SCINTILLATORCOUNTERS_H
#define SCINTILLATORCOUNTERS_H

class ScintillatorCounters
{
private:
    double location; // z-coordinate of the scintillator counter in cm (assuming it is an infinite plane in the x and y directions)
    bool direction;  // Indicates whether the detector detects x (true) or y (false) coordinates

public:
    /* BEGIN Constructor & Destructor */
    ScintillatorCounters(double location, bool direction);
    ~ScintillatorCounters();
    /* END Constructor & Destructor */

    /* BEGIN Getters */
    // returns the z-coordinate of the scintillator counter in cm (assuming it is an infinite plane in the x and y directions)
    double getLocation() const;

    // returns whether the detector detects x (true) or y (false) coordinates
    bool getDirection() const;
    /* END Getters */
};

#endif /* SCINTILLATORCOUNTERS_H */
