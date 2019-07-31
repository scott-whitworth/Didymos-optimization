#ifndef orbitalMotion_h
#define orbitalMotion_h
#include "elements.h"

// Input for trajectory and trajectoryPrint:
//  Variables to be optimized:
//      x[0]-x[8]: gamma coefficients used to calculate fourier series
//      x[9]-x[11]: tau coefficients used to calculate fourier series
//      x[12]: alpha - launch angle (declination) position 
//      x[13]: beta - launch angle (declination) velocity 
//      x[14]: trip time - total time from launch to impact, sets the initial earth position
// Only tracks the last set of optimized variables 
double trajectory( double x[]);
// Writes the array of values from optimization to a binary file
double trajectoryPrint (double x[], double & lastStep, double & cost, int j, elements<double> & yOut);


// Input for earthInitial:
// x[14]: trip time - total time from launch to impact, sets the initial earth position
// calculates the earth's initial conditions for launch date based on impact date (oct. 5, 2022) minus the optimized trip time
elements<double> earthInitial(double tripTime);



#include "orbitalMotion.cpp"
#endif