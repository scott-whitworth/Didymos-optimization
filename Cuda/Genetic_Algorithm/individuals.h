#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include "../Runge_Kutta/rkParameters.h"
#include "../Motion_Eqns/elements.h"
#include "../Earth_calculations/earthInfo.h"
#include "../Config_Constants/config.h"

// Individual is a structure member of the genetric algorithm's population and has set of parameters and the resulting position and velocity
struct Individual {
    rkParameters<double> startParams; // input parameters for the run- the unique identifiers of this individual

    elements<double> finalPos; // final position of the spacecraft at end of run

    double posDiff; // difference in position between spacecraft and center of asteroid at end of run
    double velDiff; // difference in velocity between spacecraft and asteroid at end of run
    double cost;    // cost value of the individual, something that the genetic algorithm is attempting to minimize

    // Set the inital position of the spacecraft according to this Individual's launch time
    // Input: cConstants - to access v_escape value that determines component velocities
    //        launchCon - access earth element at this individuals tripTime offset in determining position and velocity
    // Output: this individuals startParams.y0 is set to the initial position and velocity of the spacecraft
    void initialize(const cudaConstants* cConstants);

    // Calculates a cost value to quantitatively evaluate this Individual
    // Input: cConstants in accessing properties such as pos_threshold, c3energy, and v_impact
    // Output: Assigns and returns this individuals cost value
    double getCost(const cudaConstants* cConstants);

    // Comparison operators overloaded to compare individuals by their cost values (determined in getCost())
    bool operator>(Individual &other);
    bool operator<(Individual &other);
    bool operator==(Individual &other);
};

#include "individuals.cpp"

#endif
