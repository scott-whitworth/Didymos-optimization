#ifndef individuals_h
#define individuals_h

#include "individuals.h"

// An experimental equation to determine cost (currently is being minimized in the genetic algorithm)
// Currently simply returns the positional difference, but could be more elaborate by adjusting the value of cost that is returned
double Individual::getCost(cudaConstants* cConstants) {
    if (this->posDiff < cConstants->pos_threshold) {
        cost = abs(this->velDiff - cConstants->v_impact)/cConstants->c3energy;
    }
    else {
        cost = posDiff * cConstants->c3energy;
    }
    return cost;
}

bool Individual::operator>(Individual &other) {
    if (this->cost > other.cost) {
        return true;
    }
    else {
        return false;
    }
}

bool Individual::operator<(Individual &other) {
    if (this->cost < other.cost) {
        return true;
    }
    else {
        return false;
    }
}

bool Individual::operator==(Individual &other) {
    if (this->cost == other.cost) {
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
}
*/

// Initialize's the Individual's location and velocity based on earth's location/velocity at starting trip time
void Individual::initialize(const cudaConstants* cConstants) {
    elements<double> earth = launchCon->getCondition(this->startParams.tripTime); //get Earth's position and velocity at launch

    this->startParams.y0 = elements<double>( // calculate the starting position and velocity of the spacecraft from Earth's position and velocity and spacecraft launch angles
        earth.r+ESOI*cos(this->startParams.alpha),
        earth.theta+asin(sin(M_PI-this->startParams.alpha)*ESOI/earth.r),
        earth.z, // The spacecraft Individual is set to always be in-plane (no initial Z offset relative to earth) 
        earth.vr+cos(this->startParams.zeta)*sin(this->startParams.beta)*cConstants->v_escape, 
        earth.vtheta+cos(this->startParams.zeta)*cos(this->startParams.beta)*cConstants->v_escape,
        earth.vz+sin(this->startParams.zeta)*cConstants->v_escape);
}

#endif