#define _USE_MATH_DEFINES // for use of M_PI
#include "runge_kuttaCUDA.cuh"
#include "runge_kutta.h" // used for rkCalc()
#include "acceleration.h" //used for calc_accel() and calc_coast()
#include "rkParameters.h" // the struct containing the values passed to rk4simple()
#include "orbitalMotion.h"
#include <math.h>
#include <iostream>

void callRK(){
    std::cout << "top of callRK()" << std::endl;
    int numThreads = 1;

    double timeInitial = 0;
    double absTol = RK_TOL;
    double stepSize = (orbitalPeriod - timeInitial) / MAX_NUMSTEPS;

    // example values
    /*-------------------------------------------------------------------------------------*/
    double gamma[] = {10, 10, 10, 10, 10, 10, 10, 10, 10};
    double tau[] = {3, 3, 3, 3, 3};
    double coast[] = {5, 5, 5, 5, 5};

    elements<double> finalPos;

    elements<double> earth = earthInitial(2.0);
    // timeFinal, accel, wetMass, 
    // r, theta, z, 
    // vr, vtheta, vz, 
    // gamma[], tau[], coast[], coastThreshold0
    rkParameters<double> example
    (2.0, 0.0, WET_MASS, 
    earth.r+ESOI*cos(10), earth.theta+asin(sin(M_PI-10)*ESOI/earth.r), earth.z,
    earth.vr+sin(3)*vEscape, earth.vtheta+cos(3)*vEscape, earth.vz,
    gamma, tau, coast, 0.05);
    
    rkParameters<double> inputParameters[] = {example};
    /*-------------------------------------------------------------------------------------*/

    rkParameters<double> *devInputParameters; 
    double *devTimeInitial;
    double *devStepSize;
    double *devAbsTol;
    elements<double> *devFinalPos;
   

    cudaMalloc((void**) &devInputParameters, numThreads * sizeof(rkParameters<double>));
    cudaMalloc((void**) &devTimeInitial, sizeof(double));
    cudaMalloc((void**) &devStepSize, sizeof(double));
    cudaMalloc((void**) &devAbsTol, sizeof(double));
    cudaMalloc((void**) &devFinalPos, sizeof(elements<double>));

    cudaMemcpy(devInputParameters, inputParameters, numThreads * sizeof(rkParameters<double>), cudaMemcpyHostToDevice);
    cudaMemcpy(devTimeInitial, &timeInitial, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(devStepSize, &stepSize, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(devAbsTol, &absTol, sizeof(double), cudaMemcpyHostToDevice);

    rk4SimpleCUDA<<<1,1>>>(devInputParameters, devTimeInitial, devStepSize, devAbsTol, devFinalPos);

    cudaMemcpy(&finalPos, devFinalPos, sizeof(elements<double>), cudaMemcpyDeviceToHost);
    
    cudaFree(devInputParameters);
    cudaFree(devTimeInitial);
    cudaFree(devStepSize);
    cudaFree(devAbsTol);
    
    elements<double> rk4SimpleOutput;
    example.parametersRK4Simple(timeInitial, stepSize, absTol, rk4SimpleOutput);

    std::cout << finalPos << std::endl;
    std::cout << rk4SimpleOutput << std::endl;
}

// seperate conditions are passed for each thread, but timeInitial, stepSize, and absTol are the same for every thread
__global__ void rk4SimpleCUDA(rkParameters<double> * rkParametersList, double *timeInitial, double *stepSize, double *absTol, elements<double> *finalPos){
    int threadId = threadIdx.x;
    rkParameters<double> threadRKParameters = rkParametersList[threadId]; // get the parameters for this thread

    elements<double> curPos = threadRKParameters.y0; // start with the initial conditions of the spacecraft

    elements<double> k1, k2, k3, k4, k5, k6, k7; // k variables for Runge-Kutta calculation of y based off the spacecraft's final state

    double curTime = *timeInitial; // setting time equal to the start time

    thruster<double> NEXT = thruster<double>(1); // corresponds NEXT thruster to type 1 in thruster.h

    double massFuelSpent = 0; // mass of total fuel expended (kg) starts at 0

    double deltaT; // change in time for calc_accel()

    double coast; // to hold the result from calc_coast()

    while(curTime < threadRKParameters.timeFinal){
        deltaT = *stepSize;

        coast = calc_coast(threadRKParameters.coefficients, curTime, threadRKParameters.timeFinal);
        threadRKParameters.accel = calc_accel(curPos.r, curPos.z, NEXT, massFuelSpent, deltaT, coast, threadRKParameters.wetMass);

        elements<double> v; // holds output of previous value from rkCalc

        // calculate k values and get new value of y
        rkCalc(curTime, threadRKParameters.timeFinal, *stepSize, curPos, threadRKParameters.coefficients, threadRKParameters.accel, v, curPos); 

        curTime += *stepSize; // update the current time in the simulation
        *stepSize *= calc_scalingFactor(v,curPos-v,*absTol,*stepSize); // Alter the step size for the next iteration

        // The step size cannot exceed the total time divided by 2 and cannot be smaller than the total time divided by 1000
        if (*stepSize > (threadRKParameters.timeFinal - *timeInitial) / 2){
            *stepSize = (threadRKParameters.timeFinal - *timeInitial) / 2;
        }
        else if (*stepSize < ((threadRKParameters.timeFinal - *timeInitial) / 1000)){
            *stepSize = (threadRKParameters.timeFinal - *timeInitial) / 1000;
        }

        if((curTime + *stepSize) > threadRKParameters.timeFinal)
            *stepSize = (threadRKParameters.timeFinal - curTime); // shorten the last step to end exactly at time final

        // if the spacecraft is within 0.5 au of the sun, the radial position of the spacecraft artificially increases to 1000, to force that path to not be used in the optimization.
        if (curPos.r < 0.5)
        {
            curPos.r = 1000;
        }
    }

    *finalPos = curPos;
}