#include "individuals.h"

bool greaterInd(Individual first, Individual second){
    double posRatio = getPosRatio(other);
    double firstSum = first.posDiff * posRatio + first.velDiff * (1.0 - posRatio);
    double secondSum = second.posDiff * posRatio + second.velDiff * (1.0 - posRatio);
    if(firstSum < secondSum){
        return true;
    }
    else{
        return false;
    }
}

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