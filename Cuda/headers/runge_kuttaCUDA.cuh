#ifndef runge_kuttaCUDA_cuh
#define runge_kuttaCUDA_cuh
#include "motion_equations.h" // Utility functions for calc_k()
#include "rkParameters.h" // rkParameters struct

void callRK();

__global__ void rk4SimpleCUDA(rkParameters<double> *rkParametersList, double *timeInitial, double *stepSize, double *absTol, elements<double> *finalPos);

#include "runge_kuttaCUDA.cu"
#endif