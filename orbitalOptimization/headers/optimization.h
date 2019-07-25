#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H


//  Purpose:
//    Optimizes the array containing coefficients for gamma and tau fourier series and alpha and beta angles (used in initial velocity of the spacecraft).
void optimizing (double *&start, double *step);

// Runs optimizing() with a series of different starting conditions to find different solutions
void optimizeStartConditions();

// Runs optimizing() iteratively using previous results as new starting point
void iterativeOptimize();

void writeTrajectoryToFile(double *start, double & cost);
#endif