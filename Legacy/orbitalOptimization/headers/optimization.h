#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H


//  Purpose:
//    Optimizes the array containing coefficients for gamma and tau fourier series and alpha and beta angles (used in initial velocity of the spacecraft).
void optimizing (double *&start, double *step);

// Runs optimizing() with a series of different starting conditions to find different solutions
void optimizeStartConditions(int executions);

// Runs optimizing() iteratively using previous results as new starting point
void iterativeOptimize();

// This function writes the succesfull results from nealder mead into binary files
// Input the parameter array, a declared variable coast and index to name the file
void writeTrajectoryToFile(double *start, double & cost, int i);

// Initial change in variable size based on the variable start value
// Delimits the search space
// it is just so the step is consistend between the different types of testing functions
void setStep(double step[]);

// Reads a binary file with an cVector that has the same order than the code
//  and checks the convergance and writes a solution file
void checkBinaryFile(int size);

#endif