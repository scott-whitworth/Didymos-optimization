#define _USE_MATH_DEFINES // for use of M_PI
#include "runge_kuttaCUDA.cuh"
#include "runge_kutta.h" // used for rkCalc()
#include "../Thrust_Files/acceleration.h" //used for calc_accel() and calc_coast()
#include "rkParameters.h" // the struct containing the values passed to rk4simple()
#include "../Earth_calculations/orbitalMotion.h"
#include "../Genetic_Algorithm/geneticAlgorithm.h" // selectWinners()
#include "../Genetic_Algorithm/ga_crossover.h" // crossover()
#include "../Genetic_Algorithm/gaConstants.h" // SURVIVOR_COUNT
#include <math.h>
#include <iostream>
#include <fstream> // for outputing to .csv file
#include <chrono>
#include <algorithm> // sort(), shuffle()
#include <random>

// THIS FUNCTION HAS NOT BEEN COMPLETED AND IS NOT IN USE
Individual bestChange(Individual original, double timeInitial, double stepSize, double absTol){
    Individual best = original;
    Individual cur = original;

    elements<double> output;

    double parameterChange;
    
    // get the original result
    best.startParams.parametersRK4Simple(timeInitial, stepSize, absTol, output);
    best.finalPos = output;
    best.posDiff =  sqrt(pow(R_FIN_AST - output.r, 2) + pow(THETA_FIN_AST - fmod(output.theta, 2 * M_PI), 2) + pow(Z_FIN_AST - output.z, 2));
    best.velDiff =  sqrt(pow(VR_FIN_AST - output.vr, 2) + pow(VTHETA_FIN_AST - output.vtheta, 2) + pow(VZ_FIN_AST - output.vz, 2));

    // get results for each changed variable
    // gamma
    parameterChange = 0.1;
    for(int i = 0; i < 7; i++){
        cur.startParams.coeff.gamma[i] += parameterChange;
        cur.startParams.parametersRK4Simple(timeInitial, stepSize, absTol, output);
        if(!greaterInd(best, cur)){
            best = cur;
        }
        cur.startParams.coeff.gamma[i] -= parameterChange;
    }
    //tau
    parameterChange = 0.1;
    for(int i = 0; i < 3; i++){
        cur.startParams.coeff.tau[i] += parameterChange;
        cur.startParams.parametersRK4Simple(timeInitial, stepSize, absTol, output);
        if(!greaterInd(best, cur)){
            best = cur;
        }
        cur.startParams.coeff.tau[i] -= parameterChange;
    }
    //coast

    return best;
}

double optimize(const int numThreads, const int blockThreads){
    double calcPerS = 0;
    time_t timeSeed = time(0);
    std::cout << "Time seed for this run: " << timeSeed << std::endl; // note there are other mt_rands in the code that use different seeds
    std::cout << "------------------------------------------------------------------------" << std::endl;
    std::mt19937_64 mt_rand(timeSeed);

     // input parameters for rk4Simple which are the same for each thread
    double timeInitial = 0; // the starting time of the trip is always defined as zero   
    double absTol = RK_TOL; // the tolerance is a constant number that is shared amongst all runs
    double stepSize = (orbitalPeriod - timeInitial) / MAX_NUMSTEPS; // the starting step size- same for each run- note that the current step size varies throughout each run


    Individual *inputParameters = new Individual[numThreads]; // contains all input parameters besides those which are always common amongst every thread

    const int numStarts = 14; // the number of different sets of starting parameters in the input file

    std::ifstream starts;
    starts.open("../optimizedVector.bin", std::ifstream::in|std::ios::binary); // a file containing the final parameters of converged results from CPU calculations

    double startDoubles;

    // sort the data into 2 dimensions
    // one row is one set of starting parameters
    // each column is a specific variable:
    //    0-6 gamma
    //    7-9 tau
    //    10-12 launch angles
    //    13 trip time
    //    14-19 coast
    double arrayCPU[numStarts][OPTIM_VARS];
    for(int i = 0; i < OPTIM_VARS; i++){ // rows
        for(int j = 0; j < numStarts; j++){ // columns
            starts.read( reinterpret_cast<char*>( &startDoubles ), sizeof startDoubles );
            arrayCPU[j][i] = startDoubles;
        }
    }
    starts.close();

     // set every thread's input parameters to a set of final values from CPU calculations for use as a good starting point
    for(int i = 0; i < numThreads; i++){
        int row = mt_rand() % numStarts;

        double tripTime = arrayCPU[row][13];

        double alpha = arrayCPU[row][10];

        double beta = arrayCPU[row][11];

        double zeta = arrayCPU[row][12];

        coefficients<double> testcoeff;
        for(int j = 0; j < testcoeff.gammaSize; j++){
        testcoeff.gamma[j] = arrayCPU[row][j];
        }

        for(int j = 0; j < testcoeff.tauSize; j++){
        testcoeff.tau[j] =  arrayCPU[row][j+7];
        }

        for(int j = 0; j < testcoeff.coastSize; j++){
        testcoeff.coast[j] = arrayCPU[row][j+14];
        }

        rkParameters<double> example(tripTime, alpha, beta, zeta, testcoeff); 

        inputParameters[i].startParams = example;
    }

    // set every thread's input parameters to random values within a reasonable range
    /*
    for(int i = 0; i < numThreads; i++){ 
        double tripTime = 365*24*3600*(std::rand() % 10001 / 10000.0 + 1.0);
        double alpha = (mt_rand() % 629) / 100.0 - 3.14;
        double beta = (mt_rand() % 629) / 100.0 - 3.14;
        double zeta = (mt_rand() % 315) / 100.0 - 1.57;

        coefficients<double> testcoeff;
        for(int j = 0; j < testcoeff.gammaSize; j++){
            testcoeff.gamma[j] = mt_rand() % 201/10.0 - 10.0;
        }
        for(int j = 0; j < testcoeff.tauSize; j++){
            testcoeff.tau[j] = mt_rand() % 201/10.0 - 10.0;
        }
        for(int j = 0; j < testcoeff.coastSize; j++){
            testcoeff.coast[j] = mt_rand() % 201/10.0 - 10.0;
        }
    
        rkParameters<double> example(tripTime, alpha, beta, zeta, testcoeff); 

        inputParameters[i].startParams = example;
    }
    */


    Individual *survivors = new Individual[SURVIVOR_COUNT]; // stores the winners of the head-to-head competition
    int newInd = numThreads; // the whole population is new the first time through the loop

    // setup output of results
    std::ofstream individualDifference;
    individualDifference.open("individualDifference.csv");
    individualDifference << "posDiff" << "," << "velDiff" << "," << "r" << "," << "theta" << "," << "z" << "," << "vr" << "," << "vtheta" << "," << "vz" << "\n";
    
    for(int i = 0; i < generationsNum; i++){
        // initialize positions for the new individuals starting at the index of the first new one and going to the end of the array
        initializePosition(inputParameters + (numThreads - newInd), newInd);

        callRK(newInd, blockThreads, inputParameters + (numThreads - newInd), timeInitial, stepSize, absTol, calcPerS); // calculate trajectories for new individuals

        for(int k = 0; k < numThreads; k++){ // if we got bad results reset the Individual to random starting values (it may still be used for crossover) 
                                            // and set the final position to be way off so it gets replaced by a new Individual
            if(isnan(inputParameters[k].finalPos.r) || isnan(inputParameters[k].finalPos.theta) || isnan(inputParameters[k].finalPos.z) 
             || isnan(inputParameters[k].finalPos.vr) || isnan(inputParameters[k].finalPos.vtheta) || isnan(inputParameters[k].finalPos.vz)){
                
                std::cout << std::endl << std::endl << "NAN FOUND" << std::endl << std::endl;

                double tripTime = 365*24*3600*(std::rand() % 10001 / 10000.0 + 1.0);
                double alpha = (mt_rand() % 629) / 100.0 - 3.14;
                double beta = (mt_rand() % 629) / 100.0 - 3.14;
                double zeta = (mt_rand() % 315) / 100.0 - 1.57;
        
                coefficients<double> testcoeff;
                for(int j = 0; j < testcoeff.gammaSize; j++){
                    testcoeff.gamma[j] = mt_rand() % 201/10.0 - 10.0;
                }
                for(int j = 0; j < testcoeff.tauSize; j++){
                    testcoeff.tau[j] = mt_rand() % 201/10.0 - 10.0;
                }
                for(int j = 0; j < testcoeff.coastSize; j++){
                    testcoeff.coast[j] = mt_rand() % 201/10.0 - 10.0;
                }
            
                rkParameters<double> example(tripTime, alpha, beta, zeta, testcoeff); 
        
                inputParameters[k].startParams = example;

                inputParameters[k].posDiff = 1.0e10;
                inputParameters[k].velDiff = 0.0;
             }
        }

        std::shuffle(inputParameters, inputParameters + numThreads, mt_rand); // shuffle the Individiuals to use random members for the competition

        selectWinners(inputParameters, SURVIVOR_COUNT, survivors); // Choose which individuals are in survivors, not necessarrily only the best ones

        std::sort(inputParameters, inputParameters + numThreads, greaterInd); // put the individuals in order so we can replace the worst ones

        // finding the best variable to change in the best Individual
        // bestChange() TO BE USED HERE


        // Display a '.' to the terminal to show that a generation has been calculated
        // if it is not the 50th generation this serves to show that a generation was calculated and survivors selected
        // This also serves to visually seperate the generation display on the terminal screen
        std::cout << '.';

        // Display and print Individuals' pos and vel difference every 50 generations to terminal and .csv file
        if((i+1) % 50 == 0)
        {   
            // Display the best and worst Individuals in this generation
            std::cout << '\n';
            std::cout << "generation: " << i << std::endl;
            std::cout << "best:" << std::endl;
            std::cout << "\tposDiff: " << inputParameters[0].posDiff << std::endl;
            std::cout << "\tvelDiff: " << inputParameters[0].velDiff << std::endl;
            std::cout << "worst:" << std::endl;
            std::cout << "\tposDiff: " << inputParameters[numThreads - 1].posDiff << std::endl;
            std::cout << "\tvelDiff: " << inputParameters[numThreads - 1].velDiff << std::endl;

            // Append every Individual into a csv file to view progress over generations
            for(int j = 0; j < numThreads; j++)
            {
                individualDifference << inputParameters[j].posDiff << ","  << inputParameters[j].velDiff << ","
                << inputParameters[j].finalPos.r << "," << inputParameters[j].finalPos.theta << "," << inputParameters[j].finalPos.z << ","
                << inputParameters[j].finalPos.vr << "," << inputParameters[j].finalPos.vtheta << "," << inputParameters[j].finalPos.vz << "," << "\n";
            }
            individualDifference << "\n";
        }

        // the annnealing rate passed in is scaled between ANNEAL_MAX and ANNEAL_MIN depending on which generation this is
        double new_anneal =  ANNEAL_MAX - static_cast<double>(i) / (generationsNum - 1) * (ANNEAL_MAX - ANNEAL_MIN);

        newInd = crossover(survivors, inputParameters, SURVIVOR_COUNT, numThreads, new_anneal);
    }


    // output the best Individuals of the final generation, using writeTrajectoryToFile()
    // Files outputted allows plotting of solutions in matlab
    double *start = new double[OPTIM_VARS];
    double cost = 0;
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < inputParameters[i].startParams.coeff.gammaSize; j++){
            start[GAMMA_OFFSET + j] = inputParameters[i].startParams.coeff.gamma[j];
        }
        for(int j = 0; j < inputParameters[i].startParams.coeff.tauSize; j++){
            start[TAU_OFFSET + j] = inputParameters[i].startParams.coeff.tau[j];
        }
        for(int j = 0; j < inputParameters[i].startParams.coeff.coastSize; j++){
            start[COAST_OFFSET + j] = inputParameters[i].startParams.coeff.coast[j];
        }
        start[TRIPTIME_OFFSET] = inputParameters[i].startParams.tripTime;
        start[ALPHA_OFFSET] = inputParameters[i].startParams.alpha;
        start[BETA_OFFSET] = inputParameters[i].startParams.beta;
        start[ZETA_OFFSET] = inputParameters[i].startParams.zeta;

        cost = inputParameters[i].posDiff; // just look at position difference here for now
        // could instead use a ratio between position and velocity differnce as done in comparison of Individuals
        writeTrajectoryToFile(start, cost, i + 1);
    }


    individualDifference.close();

    delete [] inputParameters;
    delete [] survivors;

    return calcPerS;
}

void callRK(const int numThreads, const int blockThreads, Individual *generation, double timeInitial, double stepSize, double absTol, double & calcPerS){
    
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
    rk4SimpleCUDA<<<(numThreads+blockThreads-1)/blockThreads,blockThreads>>>(devGeneration, devTimeInitial, devStepSize, devAbsTol, numThreads);
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
__global__ void rk4SimpleCUDA(Individual *individuals, double *timeInitial, double *startStepSize, double *absTolInput, int n){
    int threadId = threadIdx.x + blockIdx.x * blockDim.x;
    if(threadId < n)
    {
        rkParameters<double> threadRKParameters = individuals[threadId].startParams; // get the parameters for this thread

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

        bool coast; // to hold the result from calc_coast()

        elements<double> error; // holds output of previous value from rkCalc

        while(curTime < threadRKParameters.tripTime){

            coast = calc_coast(threadRKParameters.coeff, curTime, threadRKParameters.tripTime);
            curAccel = calc_accel(curPos.r, curPos.z, NEXT, massFuelSpent, stepSize, coast, static_cast<double>(WET_MASS));
            //curAccel = 0.;

            // calculate k values and get new value of y
            rkCalc(curTime, threadRKParameters.tripTime, stepSize, curPos, threadRKParameters.coeff, curAccel, error, k1, k2, k3, k4, k5, k6, k7); 

            curTime += stepSize; // update the current time in the simulation
            
            stepSize *= calc_scalingFactor(curPos-error,error,absTol,stepSize); // Alter the step size for the next iteration

            // The step size cannot exceed the total time divided by 2 and cannot be smaller than the total time divided by 1000
            if (stepSize > (threadRKParameters.tripTime - startTime) / 100){
                stepSize = (threadRKParameters.tripTime - startTime) / 100;
            }
            else if (stepSize < (threadRKParameters.tripTime - startTime) / 1000){
                stepSize = (threadRKParameters.tripTime - startTime) / 1000;
            }
            
            if((curTime + stepSize) > threadRKParameters.tripTime){
                stepSize = (threadRKParameters.tripTime - curTime); // shorten the last step to end exactly at time final
            }

            // if the spacecraft is within 0.5 au of the sun, the radial position of the spacecraft artificially increases to 1000, to force that path to not be used in the optimization.
            if (sqrt(pow(curPos.r,2)+pow(curPos.z,2)) < 0.5)
            {
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

__host__ void initializePosition(Individual *individuals, int size){
    for(int i=0; i<size ;i++){
        individuals[i].initialize();
    }
}

//testing functions
//------------------------------------------------------------------------------------------------------------------------------------------------------------
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
    double tripTime = 2.0;
    double stepSize = 4500.0;
    double accel = 5.0e-16;
    elements<double> *v = new elements<double>[n];
    elements<double> k1, k2, k3, k4, k5, k6, k7;


    double *devCurTime;
    double *devtripTime;
    double *devStepSize;
    double *devAccel;
    int *devN;
    elements<double> *devV;
    elements<double> *devCurPos;
    coefficients<double> *devTestCoeff;

    cudaMalloc((void**) &devCurTime, sizeof(double));
    cudaMalloc((void**) &devtripTime, sizeof(double));
    cudaMalloc((void**) &devStepSize, sizeof(double));
    cudaMalloc((void**) &devAccel, sizeof(double));
    cudaMalloc((void**) &devN, sizeof(int));
    cudaMalloc((void**) &devV, n * sizeof(elements<double>));
    cudaMalloc((void**) &devCurPos, n * sizeof(elements<double>));
    cudaMalloc((void**) &devTestCoeff, sizeof(coefficients<double>));

    cudaMemcpy(devCurTime, &curTime, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(devN, &n, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(devtripTime, &tripTime, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(devStepSize, &stepSize, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(devAccel, &accel, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(devCurPos, curPos, n * sizeof(elements<double>), cudaMemcpyHostToDevice);
    cudaMemcpy(devTestCoeff, &testCoeff, sizeof(coefficients<double>), cudaMemcpyHostToDevice);
    std::cout << "memCpy on" << std::endl;
    rkCalcTest<<<n,1>>>(devCurTime, devtripTime, devStepSize, devTestCoeff, devAccel, devV, devCurPos, devN);
    std::cout << "kernel call" << std::endl;
    std::cout << curTime << std::endl;
    std::cout << tripTime << std::endl;
    std::cout << stepSize << std::endl;
    std::cout << accel << std::endl;
    cudaMemcpy(v, devV, n * sizeof(elements<double>), cudaMemcpyDeviceToHost);
    std::cout << "memCpy off #1" << std::endl;
    cudaMemcpy(curPos, devCurPos, n * sizeof(elements<double>), cudaMemcpyDeviceToHost); 
    std::cout << "memCpy off #2" << std::endl;

    elements<double> *hostV = new elements<double>[n];

    std::cout << curTime << std::endl;
    std::cout << tripTime << std::endl;
    std::cout << stepSize << std::endl;
    std::cout << accel << std::endl;
    //std::cout << testCoeff << std::endl;

    for(int i = 0; i < n; i++){
        std::cout << "i: " << i << std::endl;
        std::cout << hostV[i] << std::endl;
        std::cout << hostCurPos[i] << std::endl;

        rkCalc(curTime, tripTime, stepSize, hostCurPos[i], testCoeff, accel, hostV[i], k1, k2, k3, k4, k5, k6, k7);
    }

    double errorTol = 1e-10;
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
    cudaFree(devtripTime);
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

__global__ void rkCalcTest(double *curTime, double *tripTime, double *stepSize, coefficients<double> *testCoeff, double *accel, elements<double> *v, elements<double> *curPos, int *n){
    int threadId = threadIdx.x + blockIdx.x * blockDim.x;
    if(threadId < *n){
        elements<double> k1, k2, k3, k4, k5, k6, k7;
        rkCalc(*curTime, *tripTime, *stepSize, curPos[threadId], *testCoeff, *accel, v[threadId], k1, k2, k3, k4, k5, k6, k7);
    }
}