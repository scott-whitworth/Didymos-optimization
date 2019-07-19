#ifndef runge_kuttaCUDA_cuh
#define runge_kuttaCUDA_cuh
#include "motion_equations.h" // Utility functions for calc_k()
#include "rkParameters.h" // rkParameters struct

double optimize(const int numThreads, const int blockThreads);

rkParameters<double>* getNewStarts(rkParameters<double> *startParameters, elements<double> *finalPositions);

// sets up parameters and allocates memory for and then calls rk4SimpleCUDA()
elements<double>* callRK(const int numThreads, const int blockThreads, rkParameters<double> *inputParameters, double timeInitial, double stepSize, double absTol, double & calcPerS);

// the simple version of the runge_kutta algorithm, on GPU
__global__ void rk4SimpleCUDA(rkParameters<double> *rkParametersList, double *timeInitial, double *startStepSize, double *absTolInput, elements<double> *finalPos, int n);

//unit tests

__global__ void rkCalcTest(double *curTime, double *timeFinal, double *stepSize, coefficients<double> *testCoeff, double *accel, elements<double> *v, elements<double> *curPos, int *n);

#include "runge_kuttaCUDA.cu"
#endif