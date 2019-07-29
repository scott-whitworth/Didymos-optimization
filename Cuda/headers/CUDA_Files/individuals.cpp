#ifndef individuals_h
#define individuals_h

#include "individuals.h"

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
    elements<double> earth = launchCon->getCondition(this->startParams.tripTime);
    this->startParams.y0 = elements<double>(earth.r+ESOI*cos(this->startParams.alpha),earth.theta+asin(sin(M_PI-this->startParams.alpha)*ESOI/earth.r),earth.z,
  earth.vr+cos(this->startParams.zeta)*sin(this->startParams.beta)*vEscape, earth.vtheta+cos(this->startParams.zeta)*cos(this->startParams.beta)*vEscape,earth.vz+
  sin(this->startParams.zeta)*vEscape);

}

#endif
