#define _USE_MATH_DEFINES // for use of M_PI
#include "runge_kuttaCUDA.cuh"
#include "runge_kutta.h" // used for rkCalc()
#include "acceleration.h" //used for calc_accel() and calc_coast()
#include "rkParameters.h" // the struct containing the values passed to rk4simple()
#include "orbitalMotion.h"
#include <math.h>
#include <iostream>
#include <fstream> // for outputing to .csv file

double callRK(const int numThreads, const int blockThreads){

    // input parameters for rk4Simple which are the same for each thread
    double timeInitial = 0; // the starting time of the trip is always defined as zero
    double absTol = RK_TOL; // the tolerance is a constant number that is shared amongst all runs
    double stepSize = (orbitalPeriod - timeInitial) / MAX_NUMSTEPS; // the starting step size- same for each run- note that the current step size varies throughout each run

    elements<double> *finalPos = new elements<double>[numThreads]; // to store the output of final position and velocity for each run

    // reasonable example values for runge kutta algorithm
    /*-------------------------------------------------------------------------------------*/
    //double timeFinal = 2.5; // number of years the trip takes
    //double timeFinal = 2.5*365.25*24*60*60; // number of years the trip takes
    double timeFinal = 75178800-3600;
    
    
    //for setting every thread's parameters to the same values
    /*
    double gamma[] = {10, 10, 10, 10, 10, 10, 10, 10, 10};
    double tau[] = {3, 3, 3, 3, 3};
    double coast[] = {5, 5, 5, 5, 5};

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
    */

    rkParameters<double> *inputParameters = new rkParameters<double>[numThreads]; // contains all input parameters besides those which are always common amongst every thread

    for(int i = 0; i < numThreads; i++){ // set every thread's input parameters
        
        /*double gamma[] = {i, i, i, i, i, i, i, i, i};
        double tau[] = {i%4, i%4, i%4, i%4, i%4};
        double coast[] = {i%6, i%6, i%6, i%6, i%6};    
    
        elements<double> earth = earthInitial(timeFinal);
        
        inputParameters[i] = rkParameters<double>(timeFinal - i / 32, 0.0, WET_MASS,
        earth.r+ESOI*cos(i), earth.theta+asin(sin(M_PI-i)*ESOI/earth.r), earth.z,
        earth.vr+sin(i%4)*vEscape, earth.vtheta+cos(i%4)*vEscape, earth.vz,
        gamma, tau, coast, 0.005 * i);*/

        double gamma[] = {10, 10, 10, 10, 10, 10, 10, 10, 10};
        double tau[] = {5, 5, 5, 5, 5};
        double coast[] = {3, 3, 3, 3, 3};    
    
        elements<double> earth = earthInitial(timeFinal);
        
        inputParameters[i] = rkParameters<double>(timeFinal, WET_MASS,
        earth.r+ESOI*cos(0.5), earth.theta+asin(sin(M_PI-0.5)*ESOI/earth.r), earth.z,
        earth.vr+sin(0.1)*vEscape, earth.vtheta+cos(0.1)*vEscape, earth.vz,
        gamma, tau, coast, 0.05);
        
        // doesn't work
        //inputParameters[i] = example;
    }
    /*-------------------------------------------------------------------------------------*/

    
    cudaEvent_t Malloc_e, MemCpyDev_e, Kernel_e, MemCpyHost_e, MemCpyHostStop_e;
    cudaEventCreate(&Malloc_e);
    cudaEventCreate(&MemCpyDev_e);
    cudaEventCreate(&Kernel_e);
    cudaEventCreate(&MemCpyHost_e);
    cudaEventCreate(&MemCpyHostStop_e);
    

    rkParameters<double> *devInputParameters; 
    double *devTimeInitial;
    double *devStepSize;
    double *devAbsTol;
    elements<double> *devFinalPos;
   
    // allocate memory for the parameters passed to the device
    cudaEventRecord(Malloc_e);
    cudaMalloc((void**) &devInputParameters, numThreads * sizeof(rkParameters<double>));
    cudaMalloc((void**) &devTimeInitial, sizeof(double));
    cudaMalloc((void**) &devStepSize, sizeof(double));
    cudaMalloc((void**) &devAbsTol, sizeof(double));
    cudaMalloc((void**) &devFinalPos, numThreads * sizeof(elements<double>));

    // copy values of parameters passed to device onto device
    cudaEventRecord(MemCpyDev_e);
    cudaMemcpy(devInputParameters, inputParameters, numThreads * sizeof(rkParameters<double>), cudaMemcpyHostToDevice);
    cudaMemcpy(devTimeInitial, &timeInitial, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(devStepSize, &stepSize, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(devAbsTol, &absTol, sizeof(double), cudaMemcpyHostToDevice);

    // GPU version of rk4Simple()
    cudaEventRecord(Kernel_e);
    rk4SimpleCUDA<<<(numThreads+blockThreads-1)/blockThreads,blockThreads>>>(devInputParameters, devTimeInitial, devStepSize, devAbsTol, devFinalPos, numThreads);


    // copy the result of the kernel onto the host
    cudaEventRecord(MemCpyHost_e);
    cudaMemcpy(finalPos, devFinalPos, numThreads * sizeof(elements<double>), cudaMemcpyDeviceToHost);
    cudaEventRecord(MemCpyHostStop_e);
    
    // free memory from device
    cudaFree(devInputParameters);
    cudaFree(devTimeInitial);
    cudaFree(devStepSize);
    cudaFree(devAbsTol);

    // CPU version of rk4Simple()
    elements<double> *rk4SimpleOutput = new elements<double>[numThreads];
    for(int i = 0; i < numThreads; i++){
        inputParameters[i].parametersRK4Simple(timeInitial, stepSize, absTol, rk4SimpleOutput[i]);
    }

    // display final r, theta, z, vr, vtheta, and vz
    double maxError = 0.001; // how much difference is allowable between the CPU and GPU results
    bool errorFound = false;
    for(int i = 0; i < numThreads; i++){
        if(abs(rk4SimpleOutput[i].r - finalPos[i].r) > maxError){
            errorFound = true;
        }
        else if(abs(rk4SimpleOutput[i].theta - finalPos[i].theta) > maxError){
            errorFound = true;
        }
        else if(abs(rk4SimpleOutput[i].z - finalPos[i].z) > maxError){
            errorFound = true;
        }
        else if(abs(rk4SimpleOutput[i].vr - finalPos[i].vr) > maxError){
            errorFound = true;
        }
        else if(abs(rk4SimpleOutput[i].vtheta - finalPos[i].vtheta) > maxError){
            errorFound = true;
        }
        else if(abs(rk4SimpleOutput[i].vz - finalPos[i].vz) > maxError){
            errorFound = true;
        }

        if(errorFound){
            std::cout << "!!ERROR FOUND!!" << std::endl;
            std::cout << "CPU output " << i << std::endl;
            std::cout << rk4SimpleOutput[i] << std::endl;
            std::cout << "GPU output " << i << std::endl;
            std::cout << finalPos[i] << std::endl;

            errorFound = false;
        }
    }

    float mallocT, memCpyDevT, kernelT, memCpyHostT;
    
    cudaEventSynchronize(MemCpyHostStop_e);

    cudaEventElapsedTime(&mallocT, Malloc_e, MemCpyDev_e);
    cudaEventElapsedTime(&memCpyDevT, MemCpyDev_e, Kernel_e);
    cudaEventElapsedTime(&kernelT, Kernel_e, MemCpyHost_e);
    cudaEventElapsedTime(&memCpyHostT, MemCpyHost_e, MemCpyHostStop_e);
    
    double rkPerS = numThreads / (kernelT / 1000.0); // how many times the Runge Kutta algorithm ran in the kernel per second

    std::cout << "Device memory allocation time: " << mallocT << " ms" << std::endl;
    std::cout << "Device memory copy time: " << memCpyDevT << " ms" << std::endl;
    std::cout << "Host memory copy time: " << memCpyHostT << " ms" << std::endl;
    std::cout << "Kernel time: " << kernelT << " ms" << std::endl;
    std::cout << "Runge Kutta calculations per second: " << rkPerS << " /s" << std::endl;

    delete [] rk4SimpleOutput;
    delete [] finalPos;
    delete [] inputParameters;
    
    return rkPerS;
}

// seperate conditions are passed for each thread, but timeInitial, stepSize, and absTol are the same for every thread
__global__ void rk4SimpleCUDA(rkParameters<double> * rkParametersList, double *timeInitial, double *startStepSize, double *absTolInput, elements<double> *finalPos, int n){
    int threadId = threadIdx.x + blockIdx.x * blockDim.x;
    if(threadId <= n)
    {
        rkParameters<double> threadRKParameters = rkParametersList[threadId]; // get the parameters for this thread

        elements<double> curPos = threadRKParameters.y0; // start with the initial conditions of the spacecraft

        // storing copies of the input values
        double stepSize = *startStepSize;
        double absTol = *absTolInput;
        double curTime = *timeInitial;
        double startTime = *timeInitial;
        double curAccel = 0;


        elements<double> k1, k2, k3, k4, k5, k6, k7; // k variables for Runge-Kutta calculation of y based off the spacecraft's final state

        thruster<double> NEXT = thruster<double>(1); // corresponds NEXT thruster to type 1 in thruster.h

        double massFuelSpent = 0; // mass of total fuel expended (kg) starts at 0

        //double deltaT; // change in time for calc_accel()

        double coast; // to hold the result from calc_coast()

        elements<double> v; // holds output of previous value from rkCalc

        while(curTime < threadRKParameters.timeFinal){
            //deltaT = stepSize;

            coast = calc_coast(threadRKParameters.coeff, curTime, threadRKParameters.timeFinal);
            curAccel = calc_accel(curPos.r, curPos.z, NEXT, massFuelSpent, stepSize, coast, threadRKParameters.wetMass);

            // calculate k values and get new value of y
            rkCalc(curTime, threadRKParameters.timeFinal, stepSize, curPos, threadRKParameters.coeff, curAccel, v, curPos); 

            curTime += stepSize; // update the current time in the simulation
            stepSize *= calc_scalingFactor(v,curPos-v,absTol,stepSize); // Alter the step size for the next iteration

            // The step size cannot exceed the total time divided by 2 and cannot be smaller than the total time divided by 1000
            if (stepSize > (threadRKParameters.timeFinal - startTime) / 100){
                stepSize = (threadRKParameters.timeFinal - startTime) / 100;
            }
            else if (stepSize < ((threadRKParameters.timeFinal - startTime) / 10000)){
                stepSize = (threadRKParameters.timeFinal - startTime) / 10000;
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
}