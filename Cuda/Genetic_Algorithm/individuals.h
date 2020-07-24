#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include "../Runge_Kutta/rkParameters.h"
#include "../Motion_Eqns/elements.h"
#include "../Earth_calculations/earthInfo.h"
#include "../Config_Constants/config.h"

// Individual is a structure member of the genetric algorithm's population and has set of parameters and the resulting position and velocity
struct Individual {
    rkParameters<double> startParams; // input parameters for the run- the unique identifiers of this individual

    elements<double> finalPos; // final position of the spacecraft at end of run

    double posDiff; // difference in position between spacecraft and center of asteroid at end of run
    double velDiff; // difference in velocity between spacecraft and asteroid at end of run
    double cost;    // cost value of the individual, something that the genetic algorithm is attempting to minimize

    // Default constructor
    Individual();

    // Set the initial position of the spacecraft according to this Individual's launch time
    // Input: cConstants - to access v_escape value that determines component velocities
    //        launchCon - access earth element at this individuals tripTime offset in determining position and velocity
    // Output: this individual's startParams.y0 is set to the initial position and velocity of the spacecraft
    Individual(const cudaConstants* cConstants);

    // Set the initial position of the spacecraft according to the newly generated parameters
    // Input: cConstants - to access c3energy value used in getCost()
    //        newInd - struct returned by generateNewIndividual()
    // Output: this individual's startParams.y0 is set to the initial position and velocity of the spacecraft
    Individual(rkParameters<double> & newInd, const cudaConstants* cConstants);

    // Calculates a posDiff value
    // Input: cConstants in accessing properties such as r_fin_ast, theta_fin_ast, and z_fin_ast
    // Output: Assigns and returns this individual's posDiff value
    double getPosDiff(const cudaConstants* cConstants);

    // Calculates a velDiff value
    // Input: cConstants in accessing properties such as vr_fin_ast, vtheta_fin_ast, and vz_fin_ast
    // Output: Assigns and returns this individual's velDiff value
    double getVelDiff(const cudaConstants* cConstants);

    // Calculates a cost value to quantitatively evaluate this Individual
    // Input: cConstants in accessing properties such as pos_threshold, c3energy, and v_impact
    // Output: Assigns and returns this individual's cost value
    double getCost(const cudaConstants* cConstants);

    // Comparison operators overloaded to compare individuals by their cost values (determined in getCost())
    bool operator>(Individual &other);
    bool operator<(Individual &other);
    bool operator==(Individual &other);
};

#include "individuals.cpp"

#endif
