#ifndef individuals_h
#define individuals_h

#include "individuals.h"

// Set the inital position of the spacecraft according to this Individual's launch time
// Input: cConstants - to access v_escape value that determines component velocities
//        launchCon - access earth element at this individuals tripTime offset in determining position and velocity
// Output: this individuals startParams.y0 is set to the initial position and velocity of the spacecraft
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

// Calculates a cost value to quantitatively evaluate this Individual
// Input: cConstants in accessing properties such as pos_threshold, c3energy, and v_impact
// Output: Assigns and returns this individuals cost value
double Individual::getCost(const cudaConstants* cConstants) {
    if (this->posDiff < cConstants->pos_threshold) {
        this->cost = (cConstants->v_impact - this->velDiff)/cConstants->c3energy;
    }
    else {
        this->cost = this->posDiff * cConstants->c3energy;
    }
    return this->cost;
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

#endif
