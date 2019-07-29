#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include "../Runge_Kutta/rkParameters.h"
#include "../Motion_Eqns/elements.h"
#include "../constants.h"
#include "../Earth_calculations/earthInfo.h"

#define POSITION_THRESH 1.0E-8 // threshold for determining weights of position and velocity in comparisons of runs

struct Individual {
    rkParameters<double> startParams; // input parameters for the run- the unique identifiers of this individual

    elements<double> finalPos; // final position of the spacecraft at end of run

    double posDiff; // difference in position between spacecraft and center of asteroid at end of run
    double velDiff; // difference in velocity between spacecraft and asteroid at end of run

    void initialize();
};

    // utility function for greater()
    // returns the amount of weight placed on the posDiff rather than velDiff in comparison
    // output: 0.0 - 1.0
    //double getPosRatio(Individual first, Individual second);
    double getPosRatio(Individual first, Individual second);

        // compares two individuals
    // other is the Individual to be compared to
    // returns a reference to the better (more desirable output; higher fitness) Individual
    //bool greaterInd(Individual first, Individual second);
    bool greaterInd(Individual first, Individual second);

#include "individuals.cpp"

#endif