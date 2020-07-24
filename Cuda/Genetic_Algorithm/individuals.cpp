#ifndef individuals_h
#define individuals_h

#include "individuals.h"

// Default constructor
Individual::Individual() {
    this->posDiff = 1.0;
    this->velDiff = 0.0;
}

// Set the initial position of the spacecraft according to the newly generated parameters
// Input: cConstants - to access c3energy value used in getCost()
//        newInd - struct returned by generateNewIndividual()
// Output: this individual's startParams.y0 is set to the initial position and velocity of the spacecraft
Individual::Individual(rkParameters<double> & newInd, const cudaConstants* cConstants) {

    this->startParams = newInd;
    elements<double> earth = launchCon->getCondition(this->startParams.tripTime); //get Earth's position and velocity at launch

    this->startParams.y0 = elements<double>( // calculate the starting position and velocity of the spacecraft from Earth's position and velocity and spacecraft launch angles
        earth.r+ESOI*cos(this->startParams.alpha),
        earth.theta+asin(sin(M_PI-this->startParams.alpha)*ESOI/earth.r),
        earth.z, // The spacecraft Individual is set to always be in-plane (no initial Z offset relative to earth) 
        earth.vr+cos(this->startParams.zeta)*sin(this->startParams.beta)*cConstants->v_escape, 
        earth.vtheta+cos(this->startParams.zeta)*cos(this->startParams.beta)*cConstants->v_escape,
        earth.vz+sin(this->startParams.zeta)*cConstants->v_escape);
}

// Calculates a posDiff value
// Input: cConstants in accessing properties such as r_fin_ast, theta_fin_ast, and z_fin_ast
// Output: Assigns and returns this individual's posDiff value
__host__ __device__ double Individual::getPosDiff(const cudaConstants* cConstants) {
    this->posDiff = sqrt(pow(cConstants->r_fin_ast - this->finalPos.r, 2) + pow(cConstants->theta_fin_ast - fmod(this->finalPos.theta, 2 * M_PI), 2) + pow(cConstants->z_fin_ast - this->finalPos.z, 2));
    return this->posDiff;
}

// Calculates a velDiff value
// Input: cConstants in accessing properties such as vr_fin_ast, vtheta_fin_ast, and vz_fin_ast
// Output: Assigns and returns this individual's velDiff value
__host__ __device__ double Individual::getVelDiff(const cudaConstants* cConstants) {
    this->velDiff = sqrt(pow(cConstants->vr_fin_ast - this->finalPos.vr, 2) + pow(cConstants->vtheta_fin_ast - this->finalPos.vtheta, 2) + pow(cConstants->z_fin_ast - this->finalPos.vz, 2));
    return this->velDiff;
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
