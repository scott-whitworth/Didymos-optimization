#ifndef runge_kuttaCUDA_cuh
#define runge_kuttaCUDA_cuh
#include "motion_equations.h" // Utility functions for calc_k()
#include "rkParameters.h" // rkParameters struct

// sets up parameters and allocates memory for and then calls rk4SimpleCUDA()
double callRK(const int numThreads, const int blockThreads);

// the simple version of the runge_kutta algorithm, on GPU
__global__ void rk4SimpleCUDA(rkParameters<double> *rkParametersList, double *timeInitial, double *startStepSize, double *absTolInput, elements<double> *finalPos, int n);

#include "runge_kuttaCUDA.cu"
#endif