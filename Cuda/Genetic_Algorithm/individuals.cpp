#ifndef individuals_h
#define individuals_h

#include "individuals.h"

// An experimental equation to determine cost (currently is being minimized in the genetic algorithm)
// Currently simply returns the positional difference, but could be more elaborate by adjusting the value of cost that is returned
double Individual::getCost() {
    double ratio = getPosRatio();
    double cost = ratio * this->posDiff - (1.0 - ratio) * this->velDiff;
    return cost;
}

bool Individual::operator>(Individual &other) {
    if (this->getCost() > other.cost) {
        return true;
    }
    else return false;
}

bool Individual::operator<(Individual &other) {
    if (this->getCost() < other.cost) {
        return true;
    }
    else return false;
}

bool Individual::operator==(Individual &other) {
    if (this->getCost() == other.cost) {
        return true;
    }
    else return false;
}

/*
// betterInd uses < operator to compare first and second, returns true if first < second
bool betterInd(Individual first, Individual second) {
    if (first < second) {
        return true;
    }
    else {
        return false;
    }
}
*/

double Individual::getPosRatio() {
    if (this->posDiff > POS_THRESH) {
        return 1.0; // focus entirely on position because the spacecraft is very far from the asteroid
    }
    else {
        return this->posDiff / POS_THRESH; // focus more on position the greater the difference is based on linear scale
    }
}

// Initialize's the Individual's location and velocity based on earth's location/velocity at starting trip time
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