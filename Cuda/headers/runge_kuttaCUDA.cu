#define _USE_MATH_DEFINES // for use of M_PI
#include "runge_kuttaCUDA.cuh"
#include "runge_kutta.h" // used for rkCalc()
#include "acceleration.h" //used for calc_accel() and calc_coast()
#include "rkParameters.h" // the struct containing the values passed to rk4simple()
#include "orbitalMotion.h"
#include <math.h>
#include <iostream>

void callRK(){
    const int numThreads = 2048;

    // input parameters for rk4Simple which are the same for each thread
    double timeInitial = 0; // the starting time of the trip is always defined as zero
    double absTol = RK_TOL; // the tolerance is a constant number that is shared amongst all runs
    double stepSize = (orbitalPeriod - timeInitial) / MAX_NUMSTEPS; // the starting step size- same for each run- note that the current step size varies throughout each run

    elements<double> finalPos[numThreads]; // to store the output of final position and velocity for each run

    // reasonable example values for runge kutta algorithm
    /*-------------------------------------------------------------------------------------*/
    double gamma[] = {10, 10, 10, 10, 10, 10, 10, 10, 10};
    double tau[] = {3, 3, 3, 3, 3};
    double coast[] = {5, 5, 5, 5, 5};

    double timeFinal = 2.0; // number of years the trip takes

    elements<double> earth = earthInitial(timeFinal);
    // timeFinal, accel, wetMass, 
    // r, theta, z, 
    // vr, vtheta, vz, 
    // gamma[], tau[], coast[], coastThreshold0
    rkParameters<double> example // contains all input parameters besides those which are always common amongst every thread
    (timeFinal, 0.0, WET_MASS, 
    earth.r+ESOI*cos(10), earth.theta+asin(sin(M_PI-10)*ESOI/earth.r), earth.z,
    earth.vr+sin(3)*vEscape, earth.vtheta+cos(3)*vEscape, earth.vz,
    gamma, tau, coast, 0.05);
    
    rkParameters<double> inputParameters[numThreads];
    for(int i = 0; i < numThreads; i++){ // set every thread's input parameters to the same value for now
        inputParameters[i] = example;
    }
    /*-------------------------------------------------------------------------------------*/

    rkParameters<double> *devInputParameters; 
    double *devTimeInitial;
    double *devStepSize;
    double *devAbsTol;
    elements<double> *devFinalPos;
   
    // allocate memory for the parameters passed to the device
    cudaMalloc((void**) &devInputParameters, numThreads * sizeof(rkParameters<double>));
    cudaMalloc((void**) &devTimeInitial, sizeof(double));
    cudaMalloc((void**) &devStepSize, sizeof(double));
    cudaMalloc((void**) &devAbsTol, sizeof(double));
    cudaMalloc((void**) &devFinalPos, numThreads * sizeof(elements<double>));

    // copy values of parameters passed to device onto device
    cudaMemcpy(devInputParameters, inputParameters, numThreads * sizeof(rkParameters<double>), cudaMemcpyHostToDevice);
    cudaMemcpy(devTimeInitial, &timeInitial, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(devStepSize, &stepSize, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(devAbsTol, &absTol, sizeof(double), cudaMemcpyHostToDevice);

    // GPU version of rk4Simple()
    rk4SimpleCUDA<<<1024,numThreads/1024>>>(devInputParameters, devTimeInitial, devStepSize, devAbsTol, devFinalPos);

    // copy the result of the kernel onto the host
    cudaMemcpy(&finalPos, devFinalPos, numThreads * sizeof(elements<double>), cudaMemcpyDeviceToHost);
    
    // free memory from device
    cudaFree(devInputParameters);
    cudaFree(devTimeInitial);
    cudaFree(devStepSize);
    cudaFree(devAbsTol);
    
    // CPU version of rk4Simple()
    elements<double> rk4SimpleOutput;
    example.parametersRK4Simple(timeInitial, stepSize, absTol, rk4SimpleOutput);

    // display final r, theta, z, vr, vtheta, and vz
    std::cout << "CPU output:" << std::endl;
    std::cout << rk4SimpleOutput << std::endl;

    std::cout << "GPU output" << std::endl;
    for(int i = 0; i < numThreads; i+=100){
        std::cout << finalPos[i] << std::endl;
    }
}

// seperate conditions are passed for each thread, but timeInitial, stepSize, and absTol are the same for every thread
__global__ void rk4SimpleCUDA(rkParameters<double> * rkParametersList, double *timeInitial, double *startStepSize, double *absTolInput, elements<double> *finalPos){
    int threadId = threadIdx.x + blockIdx.x * blockDim.x;
    rkParameters<double> threadRKParameters = rkParametersList[threadId]; // get the parameters for this thread

    elements<double> curPos = threadRKParameters.y0; // start with the initial conditions of the spacecraft

    // storing copies of the input values
    double stepSize = *startStepSize;
    double absTol = *absTolInput;
    double curTime = *timeInitial;

    elements<double> k1, k2, k3, k4, k5, k6, k7; // k variables for Runge-Kutta calculation of y based off the spacecraft's final state

    thruster<double> NEXT = thruster<double>(1); // corresponds NEXT thruster to type 1 in thruster.h

    double massFuelSpent = 0; // mass of total fuel expended (kg) starts at 0

    double deltaT; // change in time for calc_accel()

    double coast; // to hold the result from calc_coast()

    elements<double> v; // holds output of previous value from rkCalc

    while(curTime < threadRKParameters.timeFinal){
        //deltaT = stepSize;

        coast = calc_coast(threadRKParameters.coefficients, curTime, threadRKParameters.timeFinal);
        threadRKParameters.accel = calc_accel(curPos.r, curPos.z, NEXT, massFuelSpent, stepSize, coast, threadRKParameters.wetMass);

        // calculate k values and get new value of y
        rkCalc(curTime, threadRKParameters.timeFinal, stepSize, curPos, threadRKParameters.coefficients, threadRKParameters.accel, v, curPos); 

        curTime += stepSize; // update the current time in the simulation
        stepSize *= calc_scalingFactor(v,curPos-v,absTol,stepSize); // Alter the step size for the next iteration

        // The step size cannot exceed the total time divided by 2 and cannot be smaller than the total time divided by 1000
        if (stepSize > (threadRKParameters.timeFinal - *timeInitial) / 2){
            stepSize = (threadRKParameters.timeFinal - *timeInitial) / 2;
        }
        else if (stepSize < ((threadRKParameters.timeFinal - *timeInitial) / 1000)){
            stepSize = (threadRKParameters.timeFinal - *timeInitial) / 1000;
        }

        if((curTime + stepSize) > threadRKParameters.timeFinal)
            stepSize = (threadRKParameters.timeFinal - curTime); // shorten the last step to end exactly at time final

        // if the spacecraft is within 0.5 au of the sun, the radial position of the spacecraft artificially increases to 1000, to force that path to not be used in the optimization.
        if (curPos.r < 0.5)
        {
            curPos.r = 1000;
        }
    }

    finalPos[threadId] = curPos; // output to this thread's index
}