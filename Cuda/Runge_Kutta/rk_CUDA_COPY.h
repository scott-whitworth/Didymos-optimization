#ifndef RK_CUDA_COPY_H
#define RK_CUDA_COPY_H

#include "../Runge_Kutta/runge_kutta.h"

//__global__ void rk4SimpleCUDA(Individual *individuals, double *timeInitial, double *startStepSize, double *absTolInput, int n, thruster<double> thrust, const cudaConstants* cConstant) {
void rk4sys_C_CLONE(Individual *individuals, double *timeInitial, double *startStepSize, double *absTolInput, 
                    int n, thruster<double> thrust, const cudaConstants* cConstant);


#include "rk_CUDA_COPY.cpp"

#endif