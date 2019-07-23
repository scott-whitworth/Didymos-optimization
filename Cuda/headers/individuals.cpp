#include "individual.h"
#include "constants.h"

Individual* Individual::greater(Individual *other){
    double posRatio = getPosRatio(other);
    double thisSum = this.posDiff * posRatio + this.velDiff * (1.0 - posRatio);
    double otherSum = other->posDiff * posRatio + other->velDiff * (1.0 - posRatio);
    if(thisSum < otherSum){
        return &this;
    }
    else{
        return other;
    }
}

double Individual::getPosRatio(Individual *other){
    double greaterDiff = this.posDiff; // get the greater position difference
    if(other->posDiff > greaterDiff){
        greaterDiff = other->posDiff;
    }

    if(greaterDiff > POSITION_THRESH){
        return 1.0; // focus entirely on position because the spacecraft is very far from the asteroid
    }
    else{
        return greaterDiff / POSITION_THRESH; // focus more on position the greater the difference is based on linear scale
    }
}