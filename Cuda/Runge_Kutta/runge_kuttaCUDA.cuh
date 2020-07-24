#ifndef runge_kuttaCUDA_cuh
#define runge_kuttaCUDA_cuh
#include "../Motion_Eqns/motion_equations.h" // Utility functions for calc_k()
#include "rkParameters.h" // rkParameters struct
#include "../Genetic_Algorithm/individuals.h"
#include "../Thrust_Files/thruster.h" // thrust.type
#include "../Config_Constants/config.h"

Individual* getNewStarts(Individual *prevGen);

// sets up parameters and allocates memory for and then calls rk4SimpleCUDA()
void callRK(const int numThreads, const int blockThreads, Individual *generation, double timeInitial, double stepSize, double absTol, double & calcPerS, thruster<double> thrust, const cudaConstants* cConstant);

// the simple version of the runge_kutta algorithm, on GPU
__global__ void rk4SimpleCUDA(Individual *individuals, double *timeInitial, double *startStepSize, double *absTolInput, int n, thruster<double> thrust, const cudaConstants* cConstant);

//this is used to find distinguishable difference between two positions
//Parameters:
//      p1, p2: positions that will be compared to each other
//      distinguishRate: the rate that this will divide from p1 and p2
//output: boolean true if there is is no distinguishable difference
bool distinguishableDifference(double p1, double p2, double distinguishRate);

#include "runge_kuttaCUDA.cu"
#endif