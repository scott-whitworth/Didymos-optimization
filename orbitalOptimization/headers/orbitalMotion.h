#ifndef orbitalMotion_h
#define orbitalMotion_h
#include "elements.h"

// Only tracks the last set of optimized variables 
double trajectory( double x[]);

// Writes the array of values from optimization to a binary file
double trajectoryPrint (double x[]);

elements<double> earthInitial(double tripTime);


#include "orbitalMotion.cpp"
#endif