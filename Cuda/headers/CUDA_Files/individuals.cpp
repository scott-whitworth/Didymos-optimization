#ifndef individuals_h
#define individuals_h

#include "individuals.h"
#include "gaConstants.h" // POSITION_THRESH

bool greaterInd(Individual first, Individual second){
    double posRatio = getPosRatio(first, second);
    double firstSum = first.posDiff * posRatio - first.velDiff * (1.0 - posRatio); // total cost is a mix of position and velocity
    double secondSum = second.posDiff * posRatio - second.velDiff * (1.0 - posRatio);
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

void Individual::initialize(){
    elements<double> earth = launchCon->getCondition(this->startParams.tripTime); //get Earth's position and velocity at launch

    this->startParams.y0 = elements<double>( // calculate the starting position and velocity of the spacecraft from Earth's position and velocity and spacecraft launch angles
    earth.r+ESOI*cos(this->startParams.alpha),
    earth.theta+asin(sin(M_PI-this->startParams.alpha)*ESOI/earth.r),
    earth.z,
    earth.vr+cos(this->startParams.zeta)*sin(this->startParams.beta)*vEscape, 
    earth.vtheta+cos(this->startParams.zeta)*cos(this->startParams.beta)*vEscape,
    earth.vz+sin(this->startParams.zeta)*vEscape);
}

#endif
