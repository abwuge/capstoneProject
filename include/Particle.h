#ifndef PARTICLE_H
#define PARTICLE_H

#include <iostream>
#include <vector>
#include <TVector3.h>

class Particle {

    private:
        double charge;
        double mass;
        TVector3 momentum;
        TVector3 position;

    public:

        Particle(double charge, double mass, TVector3 momentum, TVector3 position);
        ~Particle();

        // Getters
        double getCharge() const;
        double getMass() const;
        TVector3 getMomentum() const;
        TVector3 getPosition() const;

};


#endif /* PARTICLE_H */