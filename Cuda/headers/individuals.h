#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include "rkParameters.h"
#include "elements.h"
#include "constants.h"

#define POSITION_THRESH 1.0E-8 // threshold for determining weights of position and velocity in comparisons of runs

struct Individual {
    rkParameters<double> startParams; // input parameters for the run- the unique identifiers of this individual

    elements<double> finalPos; // final position of the spacecraft at end of run

    double posDiff; // difference in position between spacecraft and center of asteroid at end of run
    double velDiff; // difference in velocity between spacecraft and asteroid at end of run
};

    // utility function for greater()
    // returns the amount of weight placed on the posDiff rather than velDiff in comparison
    // output: 0.0 - 1.0
    //double getPosRatio(Individual first, Individual second);
    double getPosRatio(Individual first, Individual second){
        double greaterDiff = first.posDiff; // get the greater position difference
        if(second.posDiff > greaterDiff){
            greaterDiff = second.posDiff;
        }

        if(greaterDiff > POSITION_THRESH){
            return 1.0; // focus entirely on position because the spacecraft is very far from the asteroid
        }
        else{
            return greaterDiff / POSITION_THRESH; // focus more on position the greater the difference is based on linear scale
        }
    }

        // compares two individuals
    // other is the Individual to be compared to
    // returns a reference to the better (more desirable output; higher fitness) Individual
    //bool greaterInd(Individual first, Individual second);
    bool greaterInd(Individual first, Individual second){
        double posRatio = getPosRatio(first, second);
        double firstSum = first.posDiff * posRatio + first.velDiff * (1.0 - posRatio);
        double secondSum = second.posDiff * posRatio + second.velDiff * (1.0 - posRatio);
        if(firstSum < secondSum){
            return true;
        }
        else{
            return false;
        }
    }

#endif