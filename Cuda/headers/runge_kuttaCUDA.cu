#define _USE_MATH_DEFINES // for use of M_PI
#include "runge_kuttaCUDA.cuh"
#include "runge_kutta.h" // used for rkCalc()
#include "acceleration.h" //used for calc_accel() and calc_coast()
#include "rkParameters.h" // the struct containing the values passed to rk4simple()
#include "orbitalMotion.h"
#include <math.h>
#include <iostream>
#include <fstream> // for outputing to .csv file
#include <chrono>

double optimize(const int numThreads, const int blockThreads){
    double calcPerS = 0;

    bool maxErrorMet = false;
    elements<double> *output;

    // reasonable example values for runge kutta algorithm
    /*---------------------------------------------------------------------------------------*/
     // input parameters for rk4Simple which are the same for each thread
    double timeInitial = 0; // the starting time of the trip is always defined as zero
    //double timeFinal = 2.5*365.25*24*60*60; // number of years the trip takes
    double timeFinal = 75178800-3600;
    double absTol = RK_TOL; // the tolerance is a constant number that is shared amongst all runs
    double stepSize = (orbitalPeriod - timeInitial) / MAX_NUMSTEPS; // the starting step size- same for each run- note that the current step size varies throughout each run
    
    //for setting every thread's parameters to the same values
    coefficients<double> testcoeff;
    for(int i = 0; i < testcoeff.gammaSize; i++){
        testcoeff.gamma[i] = 10;
    }
    for(int i = 0; i < testcoeff.tauSize; i++){
        testcoeff.tau[i] = 10;
    }
    for(int i = 0; i < testcoeff.coastSize; i++){
        testcoeff.coast[i] = 10;
    }
    testcoeff.coastThreshold = .5;

    elements<double> earth = earthInitial(timeFinal);
    elements<double> spaceTest(earth.r+ESOI*cos(0.5), earth.theta+asin(sin(M_PI-0.5)*ESOI/earth.r), earth.z,
                               earth.vr+sin(0.75)*vEscape, earth.vtheta+cos(0.75)*vEscape, earth.vz);

    rkParameters<double> example(timeFinal, WET_MASS,spaceTest, testcoeff); 
    /*---------------------------------------------------------------------------------------*/


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

        /*double gamma[] = {10, 10, 10, 10, 10, 10, 10, 10, 10};
        double tau[] = {5, 5, 5, 5, 5};
        double coast[] = {3, 3, 3, 3, 3};    
    
        elements<double> earth = earthInitial(timeFinal);
        
        inputParameters[i] = rkParameters<double>(timeFinal, WET_MASS,
        earth.r+ESOI*cos(0.5), earth.theta+asin(sin(M_PI-0.5)*ESOI/earth.r), earth.z,
        earth.vr+sin(0.1)*vEscape, earth.vtheta+cos(0.1)*vEscape, earth.vz,
        gamma, tau, coast, 0.05);*/
        
        //set all inputs to the same values
        inputParameters[i] = example;
    }

    //Check to see if the input data is all the same
    for(int i = 0; i < numThreads-1; i++){
        if(!inputParameters[i].compare(inputParameters[i+1],1.0)){
            std::cout << "Things are off in the starting set" << std::endl;
        }
    }

    //while(!maxErrorMet){
    for(int i = 0; i < 1; i++){
        output = callRK(numThreads, blockThreads, inputParameters, timeInitial, stepSize, absTol, calcPerS);
        inputParameters = getNewStarts(inputParameters, output);
        delete [] output;
    }
    delete [] inputParameters;

    return calcPerS;
}

rkParameters<double>* getNewStarts(rkParameters<double> *startParameters, elements<double> *finalPositions){
    //implement genetic algorithm
    //rkParameters<double> *newParameters = new rkParameters<double>[];

    
    //delete [] startParameters;

    return startParameters;
}

elements<double>* callRK(const int numThreads, const int blockThreads, rkParameters<double> *inputParameters, double timeInitial, double stepSize, double absTol, double & calcPerS){

    elements<double> *finalPos = new elements<double>[numThreads]; // to store the output of final position and velocity for each run
    double *pDiff= new double[numThreads]; //difference in position between spacecraft and asteroid at end of trajectory
    double *vDiff = new double[numThreads]; //difference in velocity between spacecraft and asteroid at end of trajectory
    
    //events for timing functions
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
    double *devPDiff;
    double *devVDiff;


    // allocate memory for the parameters passed to the device
    cudaEventRecord(Malloc_e);
    cudaMalloc((void**) &devInputParameters, numThreads * sizeof(rkParameters<double>));
    cudaMalloc((void**) &devTimeInitial, sizeof(double));
    cudaMalloc((void**) &devStepSize, sizeof(double));
    cudaMalloc((void**) &devAbsTol, sizeof(double));
    cudaMalloc((void**) &devFinalPos, numThreads * sizeof(elements<double>));
    cudaMalloc((void**) &devPDiff, numThreads * sizeof(double));
    cudaMalloc((void**) &devVDiff, numThreads * sizeof(double));

    // copy values of parameters passed to device onto device
    cudaEventRecord(MemCpyDev_e);
    cudaMemcpy(devInputParameters, inputParameters, numThreads * sizeof(rkParameters<double>), cudaMemcpyHostToDevice);
    cudaMemcpy(devTimeInitial, &timeInitial, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(devStepSize, &stepSize, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(devAbsTol, &absTol, sizeof(double), cudaMemcpyHostToDevice);

    // GPU version of rk4Simple()
    cudaEventRecord(Kernel_e);
    std::cout << "Starting kernel with: <<<" << (numThreads+blockThreads-1)/blockThreads << "," << blockThreads << ">>>\n";
    rk4SimpleCUDA<<<(numThreads+blockThreads-1)/blockThreads,blockThreads>>>(devInputParameters, devTimeInitial, devStepSize, devAbsTol, devFinalPos, devPDiff, devVDiff, numThreads);

    // copy the result of the kernel onto the host
    cudaEventRecord(MemCpyHost_e);
    cudaMemcpy(finalPos, devFinalPos, numThreads * sizeof(elements<double>), cudaMemcpyDeviceToHost);
    cudaMemcpy(pDiff, devPDiff, numThreads * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(vDiff, devVDiff, numThreads * sizeof(double), cudaMemcpyDeviceToHost);
    cudaEventRecord(MemCpyHostStop_e);
    
    // free memory from device
    cudaFree(devInputParameters);
    cudaFree(devTimeInitial);
    cudaFree(devStepSize);
    cudaFree(devAbsTol);
    cudaFree(devPDiff);
    cudaFree(devVDiff);

    // CPU version of rk4Simple()
    // only calculate once since all input parameters are currently the same
    //elements<double> rk4SimpleOutput;
    //inputParameters[0].parametersRK4Simple(timeInitial, stepSize, absTol, rk4SimpleOutput);

    elements<double> *rk4SimpleOutput = new elements<double>[numThreads];

    auto start_timer = std::chrono::high_resolution_clock::now();

    for(int i = 0; i < numThreads; i++){
        inputParameters[i].parametersRK4Simple(timeInitial, stepSize, absTol, rk4SimpleOutput[i]);
          //std::cout << rk4SimpleOutput[i];
    }

    auto elapsed_time =  std::chrono::high_resolution_clock::now() - start_timer;
    std::cout << "CPU Calculation of " << numThreads << " RK Calculations took: " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed_time).count() << " ms" << std::endl;
    std::cout << "CPU Calculations: " << numThreads / (std::chrono::duration_cast<std::chrono::milliseconds>(elapsed_time).count()/1000.0) << " RK Calcs / second" <<  std::endl;

    // compare every GPU result with the one CPU result
    double maxError = 0.01; // how much difference is allowable between the CPU and GPU results
    bool errorFound = false;
    for(int i = 0; i < numThreads; i++){
        if(!finalPos[i].compare(rk4SimpleOutput[i],maxError)){
            errorFound = true;
        }

        if(errorFound){
            std::cout << "!!ERROR FOUND!!" << std::endl;
            std::cout << "CPU output " << i << std::endl;
            std::cout << rk4SimpleOutput[i] << std::endl;
            std::cout << "GPU output " << i << std::endl;
            std::cout << finalPos[i] << std::endl;
            std::cout << "Diff: " << std::endl;
            std::cout << finalPos[i]-rk4SimpleOutput[i] << std::endl;

            errorFound = false;
        }

        //testing
        std::cout << "final position" << finalPos[i] << std::endl;
        std::cout << "position difference" << pDiff[i] << std::endl;
        std::cout << "velocity difference" << vDiff[i] << std::endl << std::endl << std::endl;
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

    calcPerS = rkPerS;
    

    //delete [] rk4SimpleOutput;
    delete [] pDiff;
    delete [] vDiff;
    
    //return rkPerS;
    return finalPos; //make sure to delete outside of this function
}

// seperate conditions are passed for each thread, but timeInitial, stepSize, and absTol are the same for every thread
__global__ void rk4SimpleCUDA(rkParameters<double> * rkParametersList, double *timeInitial, double *startStepSize, double *absTolInput, elements<double> *finalPos, double *finalPDiff, double *finalVDiff, int n){
    int threadId = threadIdx.x + blockIdx.x * blockDim.x;
    if(threadId < n)
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
            rkCalc(curTime, threadRKParameters.timeFinal, stepSize, curPos, threadRKParameters.coeff, curAccel, v, k1, k2, k3, k4, k5, k6, k7); 

            curTime += stepSize; // update the current time in the simulation
            stepSize *= calc_scalingFactor(curPos-v,v,absTol,stepSize); // Alter the step size for the next iteration

            // The step size cannot exceed the total time divided by 2 and cannot be smaller than the total time divided by 1000
            if (stepSize > (threadRKParameters.timeFinal - startTime) / 100){
                stepSize = (threadRKParameters.timeFinal - startTime) / 100;
            }
            else if (stepSize < ((threadRKParameters.timeFinal - startTime) / 10000)){
                stepSize = (threadRKParameters.timeFinal - startTime) / 10000;
            }

            if((curTime + stepSize) > threadRKParameters.timeFinal){
                stepSize = (threadRKParameters.timeFinal - curTime); // shorten the last step to end exactly at time final
            }

            // if the spacecraft is within 0.5 au of the sun, the radial position of the spacecraft artificially increases to 1000, to force that path to not be used in the optimization.
            if (curPos.r < 0.5)
            {
                curPos.r = 1000;
            }
        }
        finalPos[threadId] = curPos; // output to this thread's index
        finalPDiff[threadId] =  sqrt(pow(R_FIN_AST - curPos.r, 2) + pow(THETA_FIN_AST - curPos.theta, 2) + pow(Z_FIN_AST - curPos.z, 2));
        finalVDiff[threadId] =  sqrt(pow(VR_FIN_AST - curPos.vr, 2) + pow(VTHETA_FIN_AST - curPos.vtheta, 2) + pow(VZ_FIN_AST - curPos.vz, 2));
        return;
    }
    return;
}



//testing functions
void rkCalcComparison(){
    int n = 5000;

    //parameter setup

    elements<double> *curPos = new elements<double>[n];
    elements<double> *hostCurPos = new elements<double>[n];
    for(int i = 0; i < n; i++){
        curPos[i].r = static_cast<double>(rand()%1001)/1000.0 + 0.5;
        curPos[i].theta = static_cast<double>(rand()%40000)/1000.0 - 20.0;
        curPos[i].z = static_cast<double>(rand()%200001)/10000000.0 - 0.001;
        curPos[i].vr = 0.000000018;
        curPos[i].vtheta = 0.00000021;
        curPos[i].vz = -0.000000002;

        hostCurPos[i] = curPos[i];
    }



    coefficients<double> testCoeff;
    for(int i = 0; i < testCoeff.gammaSize; i++){
        testCoeff.gamma[i] = 10;
    }
    for(int i = 0; i < testCoeff.tauSize; i++){
        testCoeff.tau[i] = 10;
    }
    for(int i = 0; i < testCoeff.coastSize; i++){
        testCoeff.coast[i] = 0.5;
    }
    testCoeff.coastThreshold = 0.05;

    double curTime = 1.0;
    double timeFinal = 2.0;
    double stepSize = 4500.0;
    double accel = 5.0e-16;
    elements<double> *v = new elements<double>[n];
    elements<double> k1, k2, k3, k4, k5, k6, k7;


    double *devCurTime;
    double *devTimeFinal;
    double *devStepSize;
    double *devAccel;
    int *devN;
    elements<double> *devV;
    elements<double> *devCurPos;
    coefficients<double> *devTestCoeff;

    cudaMalloc((void**) &devCurTime, sizeof(double));
    cudaMalloc((void**) &devTimeFinal, sizeof(double));
    cudaMalloc((void**) &devStepSize, sizeof(double));
    cudaMalloc((void**) &devAccel, sizeof(double));
    cudaMalloc((void**) &devN, sizeof(int));
    cudaMalloc((void**) &devV, n * sizeof(elements<double>));
    cudaMalloc((void**) &devCurPos, n * sizeof(elements<double>));
    cudaMalloc((void**) &devTestCoeff, sizeof(coefficients<double>));

    cudaMemcpy(devCurTime, &curTime, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(devN, &n, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(devTimeFinal, &timeFinal, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(devStepSize, &stepSize, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(devAccel, &accel, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(devCurPos, curPos, n * sizeof(elements<double>), cudaMemcpyHostToDevice);
    cudaMemcpy(devTestCoeff, &testCoeff, sizeof(coefficients<double>), cudaMemcpyHostToDevice);
    std::cout << "memCpy on" << std::endl;
    rkCalcTest<<<n,1>>>(devCurTime, devTimeFinal, devStepSize, devTestCoeff, devAccel, devV, devCurPos, devN);
    std::cout << "kernel call" << std::endl;
    std::cout << curTime << std::endl;
    std::cout << timeFinal << std::endl;
    std::cout << stepSize << std::endl;
    std::cout << accel << std::endl;
    cudaMemcpy(v, devV, n * sizeof(elements<double>), cudaMemcpyDeviceToHost);
    std::cout << "memCpy off #1" << std::endl;
    cudaMemcpy(curPos, devCurPos, n * sizeof(elements<double>), cudaMemcpyDeviceToHost); 
    std::cout << "memCpy off #2" << std::endl;

    elements<double> *hostV = new elements<double>[n];

    std::cout << curTime << std::endl;
    std::cout << timeFinal << std::endl;
    std::cout << stepSize << std::endl;
    std::cout << accel << std::endl;
    //std::cout << testCoeff << std::endl;

    for(int i = 0; i < n; i++){
        std::cout << "i: " << i << std::endl;
        std::cout << hostV[i] << std::endl;
        std::cout << hostCurPos[i] << std::endl;

        rkCalc(curTime, timeFinal, stepSize, hostCurPos[i], testCoeff, accel, hostV[i], k1, k2, k3, k4, k5, k6, k7);
    }

    double errorTol = 0.01;
    for(int i = 0; i < n; i++){
        if(abs(v[i].r - hostV[i].r) > errorTol){
            std::cout << "Thread: " << i << std::endl;
            std::cout << "GPU v: " << v[i] << std::endl;
            std::cout << "CPU v: " << hostV[i] << std::endl;
            std::cout << "difference: " << v[i] - hostV[i] << std::endl;
            std::cout << "GPU curPos: " << curPos[i] << std::endl;
            std::cout << "CPU curPos: " << hostCurPos[i] << std::endl;
            std::cout << "difference: " << curPos[i] - hostCurPos[i] << std::endl;
        }
    }
    std::cout << "done checking for errors" <<std::endl;

    cudaFree(devCurTime);
    cudaFree(devTimeFinal);
    cudaFree(devStepSize);
    cudaFree(devAccel);
    cudaFree(devN);
    cudaFree(devV);
    cudaFree(devCurPos);
    cudaFree(devTestCoeff);

    delete [] curPos;
    delete [] hostCurPos;
    delete [] hostV;
    delete [] v;
}

__global__ void rkCalcTest(double *curTime, double *timeFinal, double *stepSize, coefficients<double> *testCoeff, double *accel, elements<double> *v, elements<double> *curPos, int *n){
    int threadId = threadIdx.x + blockIdx.x * blockDim.x;
    if(threadId < *n){
        elements<double> k1, k2, k3, k4, k5, k6, k7;
        rkCalc(*curTime, *timeFinal, *stepSize, curPos[threadId], *testCoeff, *accel, v[threadId], k1, k2, k3, k4, k5, k6, k7);
    }
}