#ifndef individuals_h
#define individuals_h

#include "individuals.h"

// An experimental equation to determine cost (that would want to be minimized in the algorithm)
double Individual::getCost() {
    double cost = this->posDiff;
    //cost -= abs(this->velDiff);
    return cost;
}

bool Individual::operator>(Individual &other) {
    if (this->getCost() > other.getCost()){
        return true;
    }
    else {
        return false;
    }
}

bool Individual::operator<(Individual &other) {
    if (this->getCost() < other.getCost()){
        return true;
    }
    else {
        return false;
    }
}

bool Individual::operator==(Individual &other) {
    if (this->getCost() == other.getCost()){
        return true;
    }
    else {
        return false;
    }
}

// betterInd uses < operator to compare first and second, returns true if first < second
bool betterInd(Individual first, Individual second) {
    if (first < second) {
        return true;
    }
    else {
        return false;
    }
}

/*
double getPosRatio(Individual first, Individual second) {
    double greaterDiff = first.posDiff; // get the greater position difference
    if (second.posDiff > greaterDiff) {
        greaterDiff = second.posDiff;
    }

    if (greaterDiff > POSITION_THRESH) {
        return 1.0; // focus entirely on position because the spacecraft is very far from the asteroid
    }
    else {
        return greaterDiff / POSITION_THRESH; // focus more on position the greater the difference is based on linear scale
    }
}*/

// initialize sets the Individual's location and velocity
void Individual::initialize() {
    elements<double> earth = launchCon->getCondition(this->startParams.tripTime); //get Earth's position and velocity at launch

    this->startParams.y0 = elements<double>( // calculate the starting position and velocity of the spacecraft from Earth's position and velocity and spacecraft launch angles
        earth.r+ESOI*cos(this->startParams.alpha),
        earth.theta+asin(sin(M_PI-this->startParams.alpha)*ESOI/earth.r),
        earth.z, // The spacecraft Individual is set to always be in-plane (no initial Z offset relative to earth) 
        earth.vr+cos(this->startParams.zeta)*sin(this->startParams.beta)*vEscape, 
        earth.vtheta+cos(this->startParams.zeta)*cos(this->startParams.beta)*vEscape,
        earth.vz+sin(this->startParams.zeta)*vEscape);
}

#endif