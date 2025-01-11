// Coordinate system definition
/******************************\
*           ·------- x         *
*          /|                  *
*         / |                  *
*        y  |                  *
*           z                  *
\******************************/

#include <iostream>
#include "test.h"

using namespace std;

// Constants
const double c = 299792458; // speed of light in m/s

// Particle properties (muon)
double charge = 1; // Charge of the particle in e
double mass0 = 0.105658; // Rest mass of the particle in GeV/c^2

int main(int argc, char* argv[]) {

    // Particle properties
    double energy = 1; // Energy of the particle in GeV

    print_hello_world();

    return 0;
}