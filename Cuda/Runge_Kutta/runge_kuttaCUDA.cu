#define _USE_MATH_DEFINES // for use of M_PI
#include "runge_kuttaCUDA.cuh"
#include "runge_kutta.h" // used for rkCalc() and distinguishableDifference()
#include "../Thrust_Files/acceleration.h" //used for calc_accel() and calc_coast()
#include "rkParameters.h" // the struct containing the values passed to rk4simple()
#include "../Earth_calculations/orbitalMotion.h"
#include "../Genetic_Algorithm/geneticAlgorithm.h" // selectWinners()
#include "../Genetic_Algorithm/ga_crossover.h" // crossover()
#include "../Genetic_Algorithm/gaConstants.h" // SURVIVOR_COUNT
#include "../constants.h" // AU
#include <math.h>
#include <iostream>
#include <fstream> // for outputing to .csv file
#include <chrono>
#include <algorithm> // sort(), shuffle()
#include <random>

// THIS FUNCTION HAS NOT BEEN COMPLETED AND IS NOT IN USE
// Individual bestChange(Individual original, double timeInitial, double stepSize, double absTol) {
//     Individual best = original;
//     Individual cur = original;

//     elements<double> output;

//     double parameterChange;
    
//     // get the original result
//     best.startParams.parametersRK4Simple(timeInitial, stepSize, absTol, output);
//     best.finalPos = output;
//     best.posDiff =  sqrt(pow(R_FIN_AST - output.r, 2) + pow(THETA_FIN_AST - fmod(output.theta, 2 * M_PI), 2) + pow(Z_FIN_AST - output.z, 2));
//     best.velDiff =  sqrt(pow(VR_FIN_AST - output.vr, 2) + pow(VTHETA_FIN_AST - output.vtheta, 2) + pow(VZ_FIN_AST - output.vz, 2));

//     // get results for each changed variable
//     // gamma
//     parameterChange = 0.1;
//     for (int i = 0; i < 7; i++) {
//         cur.startParams.coeff.gamma[i] += parameterChange;
//         cur.startParams.parametersRK4Simple(timeInitial, stepSize, absTol, output);
//         if (!greaterInd(best, cur)) {
//             best = cur;
//         }
//         cur.startParams.coeff.gamma[i] -= parameterChange;
//     }
//     //tau
//     parameterChange = 0.1;
//     for (int i = 0; i < 3; i++) {
//         cur.startParams.coeff.tau[i] += parameterChange;
//         cur.startParams.parametersRK4Simple(timeInitial, stepSize, absTol, output);
//         if (!greaterInd(best, cur)) {
//             best = cur;
//         }
//         cur.startParams.coeff.tau[i] -= parameterChange;
//     }
//     //coast

//     return best;
// }

void callRK(const int numThreads, const int blockThreads, Individual *generation, double timeInitial, double stepSize, double absTol, double & calcPerS, thruster<double> thrust) {
    
    cudaEvent_t kernelStart, kernelEnd;
    cudaEventCreate(&kernelStart);
    cudaEventCreate(&kernelEnd);

    Individual *devGeneration; 
    double *devTimeInitial;
    double *devStepSize;
    double *devAbsTol;

    // allocate memory for the parameters passed to the device
    cudaMalloc((void**) &devGeneration, numThreads * sizeof(Individual));
    cudaMalloc((void**) &devTimeInitial, sizeof(double));
    cudaMalloc((void**) &devStepSize, sizeof(double));
    cudaMalloc((void**) &devAbsTol, sizeof(double));

    // copy values of parameters passed to device onto device
    cudaMemcpy(devGeneration, generation, numThreads * sizeof(Individual), cudaMemcpyHostToDevice);
    cudaMemcpy(devTimeInitial, &timeInitial, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(devStepSize, &stepSize, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(devAbsTol, &absTol, sizeof(double), cudaMemcpyHostToDevice);

    // GPU version of rk4Simple()
    cudaEventRecord(kernelStart);
    rk4SimpleCUDA<<<(numThreads+blockThreads-1)/blockThreads,blockThreads>>>(devGeneration, devTimeInitial, devStepSize, devAbsTol, numThreads, thrust);
    cudaEventRecord(kernelEnd);

    // copy the result of the kernel onto the host
    cudaMemcpy(generation, devGeneration, numThreads * sizeof(Individual), cudaMemcpyDeviceToHost);
    
    // free memory from device
    cudaFree(devGeneration);
    cudaFree(devTimeInitial);
    cudaFree(devStepSize);
    cudaFree(devAbsTol);

    float kernelT;
    
    cudaEventSynchronize(kernelEnd);

    cudaEventElapsedTime(&kernelT, kernelStart, kernelEnd);
    
    calcPerS = numThreads / (kernelT / 1000.0); // how many times the Runge Kutta algorithm ran in the kernel per second
}


// seperate conditions are passed for each thread, but timeInitial, stepSize, and absTol are the same for every thread
__global__ void rk4SimpleCUDA(Individual *individuals, double *timeInitial, double *startStepSize, double *absTolInput, int n, thruster<double> thrust) {
    int threadId = threadIdx.x + blockIdx.x * blockDim.x;
    if (threadId < n) {
        rkParameters<double> threadRKParameters = individuals[threadId].startParams; // get the parameters for this thread

        elements<double> curPos = threadRKParameters.y0; // start with the initial conditions of the spacecraft

        // storing copies of the input values
        double stepSize = *startStepSize;
        double absTol = *absTolInput;
        double curTime = *timeInitial;
        double startTime = *timeInitial;
        double curAccel = 0;

        elements<double> k1, k2, k3, k4, k5, k6, k7; // k variables for Runge-Kutta calculation of y based off the spacecraft's final state

        double massFuelSpent = 0; // mass of total fuel expended (kg) starts at 0

        bool coast; // to hold the result from calc_coast()

        elements<double> error; // holds output of previous value from rkCalc

        while (curTime < threadRKParameters.tripTime) {

            coast = calc_coast(threadRKParameters.coeff, curTime, threadRKParameters.tripTime, thrust);
            curAccel = calc_accel(curPos.r, curPos.z, thrust, massFuelSpent, stepSize, coast, static_cast<double>(WET_MASS));
            //curAccel = 0.;

            // calculate k values and get new value of y
            rkCalc(curTime, threadRKParameters.tripTime, stepSize, curPos, threadRKParameters.coeff, curAccel, error, k1, k2, k3, k4, k5, k6, k7, thrust); 

            curTime += stepSize; // update the current time in the simulation
            
            stepSize *= calc_scalingFactor(curPos-error,error,absTol,stepSize); // Alter the step size for the next iteration

            // The step size cannot exceed the total time divided by 2 and cannot be smaller than the total time divided by 1000
            if (stepSize > (threadRKParameters.tripTime - startTime) / 100) {
                stepSize = (threadRKParameters.tripTime - startTime) / 100;
            }
            else if (stepSize < (threadRKParameters.tripTime - startTime) / 1000) {
                stepSize = (threadRKParameters.tripTime - startTime) / 1000;
            }
            
            if ( (curTime + stepSize) > threadRKParameters.tripTime) {
                stepSize = (threadRKParameters.tripTime - curTime); // shorten the last step to end exactly at time final
            }

            // if the spacecraft is within 0.5 au of the sun, the radial position of the spacecraft artificially increases to 1000, to force that path to not be used in the optimization.
            if ( sqrt(pow(curPos.r,2)+pow(curPos.z,2)) < 0.5) {
                curPos.r = 1000;

                // output to this thread's index
                individuals[threadId].finalPos = curPos;
                individuals[threadId].posDiff = 1.0e10;
                individuals[threadId].velDiff =  0.0;

                return;
            }
        }

         // output to this thread's index
        individuals[threadId].finalPos = curPos;
        individuals[threadId].posDiff = sqrt(pow(R_FIN_AST - curPos.r, 2) + pow(THETA_FIN_AST - fmod(curPos.theta, 2 * M_PI), 2) + pow(Z_FIN_AST - curPos.z, 2));
        individuals[threadId].velDiff = sqrt(pow(VR_FIN_AST - curPos.vr, 2) + pow(VTHETA_FIN_AST - curPos.vtheta, 2) + pow(VZ_FIN_AST - curPos.vz, 2));

        return;
    }
    return;
}

__host__ void initializePosition(Individual *individuals, int size) {
    for (int i = 0; i < size ;i++) {
        individuals[i].initialize();
    }
}
    
bool distinguishableDifference(double p1, double p2, double distinguishRate) {

    double p1Num = trunc(p1/distinguishRate);
    double p2Num = trunc(p2/distinguishRate);
    
    std:: cout << "p1Num: " << p1Num << std::endl;
    std:: cout << "p2Num: " << p2Num << std::endl;

    if(p1Num == p2Num) {
        return true;
    } else {
        return false;
    }

}
void writeProgressToFile(std::ofstream fout, Individual* & pool, double gen, int thread, double prevAnneal) {
    fout.write((char*)&gen, sizeof(double));                                          // 1
    fout.write((char*)&pool[thread].posDiff, sizeof(double));                 // 2
    fout.write((char*)&pool[thread].velDiff, sizeof(double));                 // 3
    fout.write((char*)&pool[thread].finalPos.r, sizeof(double));              // 4
    fout.write((char*)&pool[thread].finalPos.theta, sizeof(double));          // 5        
    fout.write((char*)&pool[thread].finalPos.z, sizeof(double));              // 6
    fout.write((char*)&pool[thread].finalPos.vr, sizeof(double));             // 7
    fout.write((char*)&pool[thread].finalPos.vtheta, sizeof(double));         // 8
    fout.write((char*)&pool[thread].finalPos.vz, sizeof(double));             // 9
    fout.write((char*)&pool[thread].startParams.y0.r, sizeof(double));        // 10
    fout.write((char*)&pool[thread].startParams.y0.theta, sizeof(double));    // 11
    fout.write((char*)&pool[thread].startParams.y0.z, sizeof(double));        // 12
    fout.write((char*)&pool[thread].startParams.y0.vr, sizeof(double));       // 13
    fout.write((char*)&pool[thread].startParams.y0.vtheta, sizeof(double));   // 14
    fout.write((char*)&pool[thread].startParams.y0.vz, sizeof(double));       // 15
    fout.write((char*)&pool[thread].startParams.alpha, sizeof(double));       // 16
    fout.write((char*)&pool[thread].startParams.beta, sizeof(double));        // 17
    fout.write((char*)&pool[thread].startParams.zeta, sizeof(double));        // 18
    fout.write((char*)&prevAnneal, sizeof(double));                           // 19
    fout.write((char*)&pool[thread].startParams.tripTime, sizeof(double));    // 20
}

//testing functions
//------------------------------------------------------------------------------------------------------------------------------------------------------------
// void rkCalcComparison() {
//     int n = 5000;

//     //parameter setup

//     elements<double> *curPos = new elements<double>[n];
//     elements<double> *hostCurPos = new elements<double>[n];
//     for(int i = 0; i < n; i++) {
//         curPos[i].r = static_cast<double>(rand()%1001)/1000.0 + 0.5;
//         curPos[i].theta = static_cast<double>(rand()%40000)/1000.0 - 20.0;
//         curPos[i].z = static_cast<double>(rand()%200001)/10000000.0 - 0.001;
//         curPos[i].vr = 0.000000018;
//         curPos[i].vtheta = 0.00000021;
//         curPos[i].vz = -0.000000002;

//         hostCurPos[i] = curPos[i];
//     }



//     coefficients<double> testCoeff;
//     for(int i = 0; i < testCoeff.gammaSize; i++) {
//         testCoeff.gamma[i] = 10;
//     }
//     for(int i = 0; i < testCoeff.tauSize; i++) {
//         testCoeff.tau[i] = 10;
//     }
//     for(int i = 0; i < testCoeff.coastSize; i++) {
//         testCoeff.coast[i] = 0.5;
//     }
//     testCoeff.coastThreshold = 0.05;

//     double curTime = 1.0;
//     double tripTime = 2.0;
//     double stepSize = 4500.0;
//     double accel = 5.0e-16;
//     elements<double> *v = new elements<double>[n];
//     elements<double> k1, k2, k3, k4, k5, k6, k7;


//     double *devCurTime;
//     double *devtripTime;
//     double *devStepSize;
//     double *devAccel;
//     int *devN;
//     elements<double> *devV;
//     elements<double> *devCurPos;
//     coefficients<double> *devTestCoeff;

//     cudaMalloc((void**) &devCurTime, sizeof(double));
//     cudaMalloc((void**) &devtripTime, sizeof(double));
//     cudaMalloc((void**) &devStepSize, sizeof(double));
//     cudaMalloc((void**) &devAccel, sizeof(double));
//     cudaMalloc((void**) &devN, sizeof(int));
//     cudaMalloc((void**) &devV, n * sizeof(elements<double>));
//     cudaMalloc((void**) &devCurPos, n * sizeof(elements<double>));
//     cudaMalloc((void**) &devTestCoeff, sizeof(coefficients<double>));

//     cudaMemcpy(devCurTime, &curTime, sizeof(double), cudaMemcpyHostToDevice);
//     cudaMemcpy(devN, &n, sizeof(int), cudaMemcpyHostToDevice);
//     cudaMemcpy(devtripTime, &tripTime, sizeof(double), cudaMemcpyHostToDevice);
//     cudaMemcpy(devStepSize, &stepSize, sizeof(double), cudaMemcpyHostToDevice);
//     cudaMemcpy(devAccel, &accel, sizeof(double), cudaMemcpyHostToDevice);
//     cudaMemcpy(devCurPos, curPos, n * sizeof(elements<double>), cudaMemcpyHostToDevice);
//     cudaMemcpy(devTestCoeff, &testCoeff, sizeof(coefficients<double>), cudaMemcpyHostToDevice);
//     std::cout << "memCpy on" << std::endl;
//     rkCalcTest<<<n,1>>>(devCurTime, devtripTime, devStepSize, devTestCoeff, devAccel, devV, devCurPos, devN);
//     std::cout << "kernel call" << std::endl;
//     std::cout << curTime << std::endl;
//     std::cout << tripTime << std::endl;
//     std::cout << stepSize << std::endl;
//     std::cout << accel << std::endl;
//     cudaMemcpy(v, devV, n * sizeof(elements<double>), cudaMemcpyDeviceToHost);
//     std::cout << "memCpy off #1" << std::endl;
//     cudaMemcpy(curPos, devCurPos, n * sizeof(elements<double>), cudaMemcpyDeviceToHost); 
//     std::cout << "memCpy off #2" << std::endl;

//     elements<double> *hostV = new elements<double>[n];

//     std::cout << curTime << std::endl;
//     std::cout << tripTime << std::endl;
//     std::cout << stepSize << std::endl;
//     std::cout << accel << std::endl;
//     //std::cout << testCoeff << std::endl;

//     for (int i = 0; i < n; i++) {
//         std::cout << "i: " << i << std::endl;
//         std::cout << hostV[i] << std::endl;
//         std::cout << hostCurPos[i] << std::endl;

//         rkCalc(curTime, tripTime, stepSize, hostCurPos[i], testCoeff, accel, hostV[i], k1, k2, k3, k4, k5, k6, k7);
//     }

//     double errorTol = 1e-10;
//     for (int i = 0; i < n; i++) {
//         if ( abs(v[i].r - hostV[i].r) > errorTol) {
//             std::cout << "Thread: " << i << std::endl;
//             std::cout << "GPU v: " << v[i] << std::endl;
//             std::cout << "CPU v: " << hostV[i] << std::endl;
//             std::cout << "difference: " << v[i] - hostV[i] << std::endl;
//             std::cout << "GPU curPos: " << curPos[i] << std::endl;
//             std::cout << "CPU curPos: " << hostCurPos[i] << std::endl;
//             std::cout << "difference: " << curPos[i] - hostCurPos[i] << std::endl;
//         }
//     }
//     std::cout << "done checking for errors" <<std::endl;

//     cudaFree(devCurTime);
//     cudaFree(devtripTime);
//     cudaFree(devStepSize);
//     cudaFree(devAccel);
//     cudaFree(devN);
//     cudaFree(devV);
//     cudaFree(devCurPos);
//     cudaFree(devTestCoeff);

//     delete [] curPos;
//     delete [] hostCurPos;
//     delete [] hostV;
//     delete [] v;
// }

// __global__ void rkCalcTest(double *curTime, double *tripTime, double *stepSize, coefficients<double> *testCoeff, double *accel, elements<double> *v, elements<double> *curPos, int *n) {
//     int threadId = threadIdx.x + blockIdx.x * blockDim.x;
//     if (threadId < *n) {
//         elements<double> k1, k2, k3, k4, k5, k6, k7;
//         rkCalc(*curTime, *tripTime, *stepSize, curPos[threadId], *testCoeff, *accel, v[threadId], k1, k2, k3, k4, k5, k6, k7);
//     }
// }