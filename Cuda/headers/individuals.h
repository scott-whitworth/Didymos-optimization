#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include "rkParameters.h"
#include "elements.h"

struct Individual {
    rkParameters<double> startParams; // input parameters for the run- the unique identifiers of this individual

    elements<double> finalPos; // final position of the spacecraft at end of run

    double posDiff; // difference in position between spacecraft and center of asteroid at end of run
    double velDiff; // difference in velocity between spacecraft and asteroid at end of run

    // compares two individuals
    // other is the Individual to be compared to
    // returns a reference to the better (more desirable output; higher fitness) Individual
    Individual* greater(Individual *other);

    // utility function for greater()
    // returns the amount of weight placed on the posDiff rather than velDiff in comparison
    // output: 0.0 - 1.0
    double getPosRatio(Individual *other);
};

#endif