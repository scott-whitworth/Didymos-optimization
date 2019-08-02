#ifndef runge_kuttaCUDA_cuh
#define runge_kuttaCUDA_cuh
#include "../Motion_Eqns/motion_equations.h" // Utility functions for calc_k()
#include "rkParameters.h" // rkParameters struct
#include "../CUDA_Files/individuals.h"

double optimize(const int numThreads, const int blockThreads);

Individual* getNewStarts(Individual *prevGen);

// sets up parameters and allocates memory for and then calls rk4SimpleCUDA()
void callRK(const int numThreads, const int blockThreads, Individual *generation, double timeInitial, double stepSize, double absTol, double & calcPerS);

// the simple version of the runge_kutta algorithm, on GPU
__global__ void rk4SimpleCUDA(Individual *individuals, double *timeInitial, double *startStepSize, double *absTolInput, int n, int *testVal);

//unit tests

__global__ void rkCalcTest(double *curTime, double *timeFinal, double *stepSize, coefficients<double> *testCoeff, double *accel, elements<double> *v, elements<double> *curPos, int *n);

__host__ void initializePosition(Individual *individuals, int size);

#include "runge_kuttaCUDA.cu"
#endif