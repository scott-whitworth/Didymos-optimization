#define _USE_MATH_DEFINES // for use of M_PI
#include "runge_kuttaCUDA.cuh"
#include "runge_kutta.h" // used for rkCalc()
#include "../Thrust_Files/acceleration.h" //used for calc_accel() and calc_coast()
#include "rkParameters.h" // the struct containing the values passed to rk4simple()
#include "../Earth_calculations/orbitalMotion.h"
#include "../CUDA_Files/geneticAlgorithm.h" // selectWinners()
#include "../CUDA_Files/ga_crossover.h" // crossover()
#include "../CUDA_Files/gaConstants.h" // SURVIVOR_COUNT
#include <math.h>
#include <iostream>
#include <fstream> // for outputing to .csv file
#include <chrono>
#include <algorithm> // sort(), shuffle()
#include <random>

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
    parameterChange = 0.1;
    for(int i = 0; i < 6; i++){
        cur.startParams.coeff.gamma[i] += parameterChange;
        cur.startParams.parametersRK4Simple(timeInitial, stepSize, absTol, output);
        if(!greaterInd(best, cur)){
            best = cur;
        }
        cur.startParams.coeff.gamma[i] -= parameterChange;
    }

    return best;
}

double optimize(const int numThreads, const int blockThreads){
    double calcPerS = 0;
    //std::mt19937_64 mt_rand(time(0));
    std::mt19937_64 mt_rand(0);
    //bool maxErrorMet = false;

    // reasonable example values for runge kutta algorithm
    /*---------------------------------------------------------------------------------------*/
     // input parameters for rk4Simple which are the same for each thread
    double timeInitial = 0; // the starting time of the trip is always defined as zero
    //double tripTime = 2.5*365.25*24*60*60; // number of years the trip takes
    
    double absTol = RK_TOL; // the tolerance is a constant number that is shared amongst all runs
    double stepSize = (orbitalPeriod - timeInitial) / MAX_NUMSTEPS; // the starting step size- same for each run- note that the current step size varies throughout each run
    
    //for setting every thread's parameters to the same values
  
    /*---------------------------------------------------------------------------------------*/


    Individual *inputParameters = new Individual[numThreads]; // contains all input parameters besides those which are always common amongst every thread

    double arrayCPU[16][19] =  {{-0.11618,-0.93813,0.16112,0.66778,-0.78756,0.99569,0.4093,0.91157,-0.2663,0.60454,-0.53152,1.8883,-0.87309,3.392e+07,0.67193,0.45975,0.96308,0.35944,0.16128},
                                {1.6306,1.3656,-1.0624,0.067411,1.3817,-0.47046,-1.7251,0.76606,1.1041,-0.48066,2.9519,-1.3519,0.73574,5.5207e+07,-0.35152,-0.4905,-0.74388,-0.28909,-0.2466},
                                {-0.7665,-0.74499,0.78058,-0.57389,-0.067605,1.1937,1.529,0.7803,1.2685,-0.32151,-1.8358,-1.5837,0.89957,2.6405e+07,0.17785,0.25197,0.80231,-1.1902,-0.76714},
                                {0.8741,1.9769,1.8014,-0.76438,-0.0086669,0.67533,0.12452,-0.019435,0.74589,-0.27469,-0.74055,2.1006,-1.1899,4.1777e+07,-0.20959,-0.55904,0.21765,0.31739,0.10648},
                                {0.16945,-0.16863,-0.049979,-0.1288,0.041762,0.22557,-0.11758,0.059362,-0.66864,-0.51101,1.1111,3.7917,0.56411,5.0875e+07,-1.0008,-0.0013711,-1.0439,-0.519,-0.69473},
                                {1.3611,-0.31388,1.0975,0.67393,-1.066,0.76805,1.0808,-0.63892,0.023753,-0.49204,-0.21729,-1.1225,1.8982,2.7365e+07,0.70863,-0.54919,-0.22024,0.043853,-0.17199},
                                {0.26668,0.303,-0.24341,0.75273,0.80444,0.8867,-0.12381,-0.40879,-0.75649,-0.051487,-2.5249,1.6057,-0.61328,6.6164e+07,0.038794,0.80025,-0.1774,0.25521,-0.79408},
                                {-0.33794,-0.10588,0.11451,-0.11982,0.11674,0.020274,-0.18984,-0.22741,-0.24874,0.19079,2.0114,-2.0955,-0.68952,4.3259e+07,-0.76951,-0.080516,0.78251,-0.66553,-0.556},
                                {0.24924,-0.53743,0.03409,-0.43815,0.43498,0.18907,0.087696,0.27394,-0.31364,-0.048284,-2.1545,2.1532,0.96415,5.251e+07,0.0064018,0.22662,-0.66381,-0.49153,-0.52606},
                                {2.6971,0.33445,1.4728,-0.20598,-0.59999,0.13327,-0.28345,0.78588,-0.71313,-0.37493,-0.50709,-1.5517,1.0289,2.5667e+07,-0.30889,-0.18324,-0.83499,1.0481,0.11909},
                                {0.20526,-0.15981,0.16457,0.60735,0.32903,0.30262,0.20892,-0.1483,-0.24254,-0.27927,-0.009502,3.3074,0.65589,5.3982e+07,0.76371,-0.33244,0.57125,0.068041,0.69824},
                                {0.47673,-0.1147,0.22333,0.080363,0.64068,0.21915,0.18892,0.22581,0.34533,0.59555,0.49483,-0.39331,-1.8962,4.2546e+07,-0.20768,-0.77094,-0.37879,-0.56058,0.39925},
                                {-0.12378,1.2322,-0.12036,0.61289,-0.2552,0.37191,0.23688,0.072462,0.28729,0.50243,2.5028,1.7963,0.68139,5.5981e+07,0.021783,-0.12746,-0.60296,-0.13233,-0.53686},
                                {0.12762,-0.19833,0.44299,-0.20322,0.13699,-0.085633,0.49096,0.55722,-0.56092,0.67755,-0.45509,-2.0707,-0.62307,4.1228e+07,1.0357,0.36585,-0.77794,0.54923,0.93459},
                                {0.14021,-0.12283,0.13903,0.13796,-0.45269,-0.44853,0.20326,-0.31895,0.26559,0.94838,1.116,-2.2627,-0.3194,3.5226e+07,-0.98734,0.8642,-0.2093,1.1429,0.10095},
                                {-0.084377,0.63781,-0.084106,0.56619,-0.84042,-0.0043268,1.1472,-0.14925,-0.33214,0.25764,-2.1727,1.6719,1.273,5.1044e+07,-0.399,0.55074,0.41427,0.62746,-0.52202}};

for(int i = 0; i < numThreads; i++){ // set every thread's input parameters
    //double tripTime = 365*24*3600*(std::rand() % 10001 / 10000.0 + 1.0);
    int row = mt_rand() % 16;
    double tripTime = arrayCPU[row][13];
    //std::cout << "tripTime: " << tripTime << std::endl;
    //double alpha = (mt_rand() % 629) / 100.0 - 3.14;
    double alpha = arrayCPU[row][10];
    //std::cout << "alpha: " << alpha << std::endl;
    //double beta = (mt_rand() % 629) / 100.0 - 3.14;
    double beta = arrayCPU[row][11];
    //std::cout << "beta: " << beta << std::endl;
    //double zeta = (mt_rand() % 315) / 100.0 - 1.57;
    double zeta = arrayCPU[row][12];
    //std::cout << "zeta: " << zeta << std::endl;
    /*
    coefficients<double> testcoeff;
    for(int j = 0; j < testcoeff.gammaSize; j++){
    testcoeff.gamma[j] = mt_rand() % 201/10.0 - 10.0;
    }
    */
    coefficients<double> testcoeff;
    for(int j = 0; j < testcoeff.gammaSize; j++){
    testcoeff.gamma[j] = arrayCPU[row][j];
    //std::cout << "gamma: " << testcoeff.gamma[j] << std::endl;
    }
    /*
    for(int j = 0; j < testcoeff.tauSize; j++){
    testcoeff.tau[j] = mt_rand() % 201/10.0 - 10.0;
    }
    */
    for(int j = 0; j < testcoeff.tauSize; j++){
    testcoeff.tau[j] =  arrayCPU[row][j+7];
    //std::cout << "tau: " << testcoeff.tau[j] << std::endl;
    }

    /*
    for(int j = 0; j < testcoeff.coastSize; j++){
    testcoeff.coast[j] = mt_rand() % 201/10.0 - 10.0;
    }
    */
    for(int j = 0; j < testcoeff.coastSize; j++){
    testcoeff.coast[j] = arrayCPU[row][j+14];
    //std::cout << "coast: " << testcoeff.coast[j] << std::endl;
    }


    //elements<double> earth = earthInitial(tripTime);
    //elements<double> spaceTest(earth.r+ESOI*cos(alpha), earth.theta+asin(sin(M_PI-alpha)*ESOI/earth.r), earth.z,
    //earth.vr+cos(zeta)*sin(beta)*vEscape, earth.vtheta+cos(zeta)*cos(beta)*vEscape, earth.vz+sin(zeta)*vEscape);

    rkParameters<double> example(tripTime, alpha, beta, zeta, testcoeff); 

    inputParameters[i].startParams = example;

    //std::cout << "example: " << example << std::endl;
}

    /*
    for(int i = 0; i < numThreads; i++){ // set every thread's input parameters
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
    
        //elements<double> earth = earthInitial(tripTime);
        //elements<double> spaceTest(earth.r+ESOI*cos(alpha), earth.theta+asin(sin(M_PI-alpha)*ESOI/earth.r), earth.z,
            //earth.vr+cos(zeta)*sin(beta)*vEscape, earth.vtheta+cos(zeta)*cos(beta)*vEscape, earth.vz+sin(zeta)*vEscape);
    
        rkParameters<double> example(tripTime, alpha, beta, zeta, testcoeff); 

        inputParameters[i].startParams = example;
    }
    */


    Individual *survivors = new Individual[SURVIVOR_COUNT];
    int newInd = numThreads; // the whole population is new the first time through the loop

    // printing individual pos and vel difference data to a csv
    std::ofstream individualDifference;
    individualDifference.open("individualDifference.csv");
    individualDifference << "posDiff" << "," << "velDiff" << "," << "r" << "," << "theta" << "," << "z" << "," << "vr" << "," << "vtheta" << "," << "vz" << "\n";

    double annealMax = ANNEAL_MAX;
    double annealMin = ANNEAL_MIN;
    //while(!maxErrorMet){
    for(int i = 0; i < generationsNum; i++){
        auto start = std::chrono::high_resolution_clock::now();
        initializePosition(inputParameters + (numThreads - newInd), newInd); // initialize positions for new individuals
        
        // testing
        /*
        std::cout << "generation " << i << std::endl << std::endl;
        for(int j = 0; j < numThreads; j++){ 
            std::cout << inputParameters[j].startParams.y0 << std::endl;
        }
        */
        auto RK = std::chrono::high_resolution_clock::now();
        callRK(newInd, blockThreads, inputParameters + (numThreads - newInd), timeInitial, stepSize, absTol, calcPerS); // calculate trajectories for new individuals
        
        auto NaNCheck = std::chrono::high_resolution_clock::now();
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
        auto shuffleT = std::chrono::high_resolution_clock::now();
        std::shuffle(inputParameters, inputParameters + numThreads, mt_rand);

        auto competition = std::chrono::high_resolution_clock::now();
        selectWinners(inputParameters, SURVIVOR_COUNT, survivors);

        auto sort = std::chrono::high_resolution_clock::now();
        std::sort(inputParameters, inputParameters + numThreads, greaterInd);

        // finding the best variable to cahnge in the best Individual








        
        auto display = std::chrono::high_resolution_clock::now();
        std::cout << "generation: " << i << std::endl;
        std::cout << "best:" << std::endl;
        std::cout << "posDiff: " << inputParameters[0].posDiff << std::endl;
        std::cout << "velDiff: " << inputParameters[0].velDiff << std::endl;
        std::cout << "finalPos: " <<inputParameters[0].finalPos << std::endl;
        std::cout << "worst:" << std::endl;
        std::cout << "posDiff: " << inputParameters[numThreads - 1].posDiff << std::endl;
        std::cout << "velDiff: " << inputParameters[numThreads - 1].velDiff << std::endl;
        std::cout << "finalPos: " <<inputParameters[numThreads - 1].finalPos << std::endl << std::endl;


        // For csv file "individualDifference.csv"
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if(i % 50 == 0)
        {   
            for(int j = 0; j < numThreads; j++)
            {
                individualDifference << inputParameters[j].posDiff << ","  << inputParameters[j].velDiff << "," << inputParameters[j].finalPos.r << "," <<
                 inputParameters[j].finalPos.theta << "," << inputParameters[j].finalPos.z << "," << inputParameters[j].finalPos.vr << "," << 
                 inputParameters[j].finalPos.vtheta << "," << inputParameters[j].finalPos.vz << "," << "\n";
            }
            individualDifference << "\n";
        }
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


        auto crossoverT = std::chrono::high_resolution_clock::now();
        newInd = crossover(survivors, inputParameters, SURVIVOR_COUNT, numThreads, annealMax - static_cast<double>(i) / (generationsNum - 1) * (annealMax - annealMin));
        auto end = std::chrono::high_resolution_clock::now();


        // display timing metrics
       /* 
        std::chrono::duration<double> elapsedTime = RK - start;
        std::cout << "Execution speeds (seconds):" << std::endl;
        std::cout << "initializePosition(): " << elapsedTime.count() << std::endl;
        elapsedTime = NaNCheck - RK;
        std::cout << "callRK(): " << elapsedTime.count() << std::endl;
        elapsedTime = shuffleT - NaNCheck;
        std::cout << "NaN check: " << elapsedTime.count() << std::endl;
        elapsedTime = competition - shuffleT;
        std::cout << "shuffle(): " << elapsedTime.count() << std::endl;
        elapsedTime = sort - competition;
        std::cout << "selectWinners(): " << elapsedTime.count() << std::endl;
        elapsedTime = display - sort;
        std::cout << "sort(): " << elapsedTime.count() << std::endl;
        elapsedTime = crossoverT - display;
        std::cout << "display(): " << elapsedTime.count() << std::endl;
        elapsedTime = end - crossoverT;
        std::cout << "crossover(): " << elapsedTime.count() << std::endl << std::endl;
        */
        
    }

    individualDifference.close();

    delete [] inputParameters;
    delete [] survivors;

    return calcPerS;
}

void callRK(const int numThreads, const int blockThreads, Individual *generation, double timeInitial, double stepSize, double absTol, double & calcPerS){
    

    auto start2 = std::chrono::high_resolution_clock::now();
    //events for timing functions
    cudaEvent_t Malloc_e, MemCpyDev_e, Kernel_e, MemCpyHost_e, MemCpyHostStop_e;
    cudaEventCreate(&Malloc_e);
    cudaEventCreate(&MemCpyDev_e);
    cudaEventCreate(&Kernel_e);
    cudaEventCreate(&MemCpyHost_e);
    cudaEventCreate(&MemCpyHostStop_e);

    auto indiv = std::chrono::high_resolution_clock::now();
    Individual *devGeneration; 
    double *devTimeInitial;
    double *devStepSize;
    double *devAbsTol;
    
    
    int *devTestVal;
    int *testVal = new int[numThreads];


    auto allocating = std::chrono::high_resolution_clock::now();
    // allocate memory for the parameters passed to the device
    cudaEventRecord(Malloc_e);
    cudaMalloc((void**) &devGeneration, numThreads * sizeof(Individual));
    cudaMalloc((void**) &devTimeInitial, sizeof(double));
    cudaMalloc((void**) &devStepSize, sizeof(double));
    cudaMalloc((void**) &devAbsTol, sizeof(double));
    cudaMalloc((void**) &devTestVal, numThreads * sizeof(int));

    auto copyParam = std::chrono::high_resolution_clock::now();
    // copy values of parameters passed to device onto device
    cudaEventRecord(MemCpyDev_e);
    cudaMemcpy(devGeneration, generation, numThreads * sizeof(Individual), cudaMemcpyHostToDevice);
    cudaMemcpy(devTimeInitial, &timeInitial, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(devStepSize, &stepSize, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(devAbsTol, &absTol, sizeof(double), cudaMemcpyHostToDevice);


    auto rkSIM = std::chrono::high_resolution_clock::now();
    // GPU version of rk4Simple()
    cudaEventRecord(Kernel_e);
    //std::cout << "Starting kernel with: <<<" << (numThreads+blockThreads-1)/blockThreads << "," << blockThreads << ">>>\n";
    rk4SimpleCUDA<<<(numThreads+blockThreads-1)/blockThreads,blockThreads>>>(devGeneration, devTimeInitial, devStepSize, devAbsTol, numThreads, devTestVal);


    auto copyRes = std::chrono::high_resolution_clock::now();
    // copy the result of the kernel onto the host
    cudaEventRecord(MemCpyHost_e);
    cudaMemcpy(generation, devGeneration, numThreads * sizeof(Individual), cudaMemcpyDeviceToHost);
    cudaMemcpy(testVal, devTestVal, numThreads * sizeof(int), cudaMemcpyDeviceToHost);
    cudaEventRecord(MemCpyHostStop_e);
    

    auto freeMem = std::chrono::high_resolution_clock::now();
    // free memory from device
    cudaFree(devGeneration);
    cudaFree(devTimeInitial);
    cudaFree(devStepSize);
    cudaFree(devAbsTol);
    cudaFree(devTestVal);



    auto rkSIM_CPU = std::chrono::high_resolution_clock::now();
    // CPU version of rk4Simple()
    // only calculate once since all input parameters are currently the same
    //elements<double> rk4SimpleOutput;
    //inputParameters[0].parametersRK4Simple(timeInitial, stepSize, absTol, rk4SimpleOutput);

    /*elements<double> *rk4SimpleOutput = new elements<double>[numThreads];
    
    for(int i = 0; i < numThreads; i++){
        generation[i].startParams.parametersRK4Simple(timeInitial, stepSize, absTol, rk4SimpleOutput[i]);
          //std::cout << rk4SimpleOutput[i];
    }
    */
    /*
    std::cout << "CPU Calculation of " << numThreads << " RK Calculations took: " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed_time).count() << " ms" << std::endl;
    std::cout << "CPU Calculations: " << numThreads / (std::chrono::duration_cast<std::chrono::milliseconds>(elapsed_time).count()/1000.0) << " RK Calcs / second" <<  std::endl;
    */

   /*
    auto errorFinding = std::chrono::high_resolution_clock::now();
    // compare every GPU result with the one CPU result
    double maxError = 1e-12; // how much difference is allowable between the CPU and GPU results
    bool errorFound = false;
    for(int i = 0; i < numThreads; i++){
        if(!generation[i].finalPos.compare(rk4SimpleOutput[i],maxError)){
            errorFound = true;
        }

        if(errorFound){
            std::cout << "!!ERROR FOUND!!" << std::endl;
            std::cout << "CPU output " << i << std::endl;
            std::cout << rk4SimpleOutput[i] << std::endl;
            std::cout << "GPU output " << i << std::endl;
            std::cout << generation[i].finalPos << std::endl;
            std::cout << "Diff: " << std::endl;
            std::cout << generation[i].finalPos-rk4SimpleOutput[i] << std::endl;

            errorFound = false;
        }

        //testing
     
        std::cout << "final position" << generation[i].finalPos << std::endl;
        std::cout << "position difference" << generation[i].posDiff << std::endl;
        std::cout << "velocity difference" << generation[i].velDiff << std::endl << std::endl << std::endl;
   
    }
     */

    
    auto resultCheck = std::chrono::high_resolution_clock::now();
    float mallocT, memCpyDevT, kernelT, memCpyHostT;
    
    cudaEventSynchronize(MemCpyHostStop_e);

    cudaEventElapsedTime(&mallocT, Malloc_e, MemCpyDev_e);
    cudaEventElapsedTime(&memCpyDevT, MemCpyDev_e, Kernel_e);
    cudaEventElapsedTime(&kernelT, Kernel_e, MemCpyHost_e);
    cudaEventElapsedTime(&memCpyHostT, MemCpyHost_e, MemCpyHostStop_e);
    
    double rkPerS = numThreads / (kernelT / 1000.0); // how many times the Runge Kutta algorithm ran in the kernel per second

    /*
    std::cout << "Device memory allocation time: " << mallocT << " ms" << std::endl;
    std::cout << "Device memory copy time: " << memCpyDevT << " ms" << std::endl;
    std::cout << "Host memory copy time: " << memCpyHostT << " ms" << std::endl;
    std::cout << "Kernel time: " << kernelT << " ms" << std::endl;
    std::cout << "Runge Kutta calculations per second: " << rkPerS << " /s" << std::endl;
    */

    calcPerS = rkPerS;
    auto end = std::chrono::high_resolution_clock::now();

    
    //delete [] rk4SimpleOutput;


    // display timing metrics
    /*
    std::chrono::duration<double> elapsedTime = indiv - start2;
    std::cout << "Execution speeds (seconds):" << std::endl;
    std::cout << "start: " << elapsedTime.count() << std::endl;

    elapsedTime = allocating - indiv;
    std::cout << "individuals: " << elapsedTime.count() << std::endl;

    elapsedTime = copyParam - allocating;
    std::cout << "allocating: " << elapsedTime.count() << std::endl;

    elapsedTime = rkSIM - copyParam;
    std::cout << "Copy parameters: " << elapsedTime.count() << std::endl;

    elapsedTime = copyRes - rkSIM;
    std::cout << "RK simple: " << elapsedTime.count() << std::endl;

    elapsedTime = freeMem - copyRes;
    std::cout << "copy results: " << elapsedTime.count() << std::endl;

    elapsedTime = rkSIM_CPU - freeMem;
    std::cout << "Freeing memory: " << elapsedTime.count() << std::endl;

    elapsedTime = resultCheck - rkSIM_CPU;
    std::cout << "RK simple CPU: " << elapsedTime.count() << std::endl;

    elapsedTime = end - resultCheck;
    std::cout << "result check: " << elapsedTime.count() << std::endl << std::endl;
    */

    //for(int i = 0; i < numThreads; i++){
    //for(int i = 0; i < 200; i++){
    //   std::cout << "iterations before NAN: " << testVal[i] << "\terror.vr: " << generation[i].finalPos.r << "\tprevious error.vr: " << generation[i].finalPos.theta << std::endl;
    //}

    delete [] testVal;
}



// seperate conditions are passed for each thread, but timeInitial, stepSize, and absTol are the same for every thread
__global__ void rk4SimpleCUDA(Individual *individuals, double *timeInitial, double *startStepSize, double *absTolInput, int n, int *testVal){
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

        int iteration = 0;
        bool stop = false;


        elements<double> k1, k2, k3, k4, k5, k6, k7; // k variables for Runge-Kutta calculation of y based off the spacecraft's final state

        thruster<double> NEXT = thruster<double>(1); // corresponds NEXT thruster to type 1 in thruster.h

        double massFuelSpent = 0; // mass of total fuel expended (kg) starts at 0

        //double deltaT; // change in time for calc_accel()

        bool coast; // to hold the result from calc_coast()

        elements<double> error; // holds output of previous value from rkCalc
        elements<double> testEl2 = curPos;

        double startStep = stepSize;

        while(curTime < threadRKParameters.tripTime && !stop){
            //deltaT = stepSize;

            coast = calc_coast(threadRKParameters.coeff, curTime, threadRKParameters.tripTime);
            curAccel = calc_accel(curPos.r, curPos.z, NEXT, massFuelSpent, stepSize, coast, static_cast<double>(WET_MASS));
            //curAccel = 0.;

            // calculate k values and get new value of y
            rkCalc(curTime, threadRKParameters.tripTime, stepSize, curPos, threadRKParameters.coeff, curAccel, error, k1, k2, k3, k4, k5, k6, k7); 

            curTime += stepSize; // update the current time in the simulation
            
            stepSize *= calc_scalingFactor(curPos-error,error,absTol,stepSize)/2; // Alter the step size for the next iteration

            // The step size cannot exceed the total time divided by 2 and cannot be smaller than the total time divided by 1000
            if (stepSize > (threadRKParameters.tripTime - startTime) / 1000){
                stepSize = (threadRKParameters.tripTime - startTime) / 1000;
            }
            else if (stepSize < (threadRKParameters.tripTime - startTime) / 10000){
                stepSize = (threadRKParameters.tripTime - startTime) / 10000;
            }
            
            if((curTime + stepSize) > threadRKParameters.tripTime){
                stepSize = (threadRKParameters.tripTime - curTime); // shorten the last step to end exactly at time final
            }

            // if the spacecraft is within 0.5 au of the sun, the radial position of the spacecraft artificially increases to 1000, to force that path to not be used in the optimization.
            if (curPos.r < 0.5)
            {
                curPos.r = 1000;
            }
            
            
            if(isnan(curPos.r) || isnan(curPos.theta) || isnan(curPos.z) 
             || isnan(curPos.vr) || isnan(curPos.vtheta) || isnan(curPos.vz)){
                 stop = true;
             }

            if(!stop){
                iteration++;
            }
        }
        individuals[threadId].finalPos = curPos; // output to this thread's index
        testVal[threadId] = iteration;

        individuals[threadId].posDiff = pow(R_FIN_AST - curPos.r, 2) + pow(THETA_FIN_AST - fmod(curPos.theta, 2 * M_PI), 2) + pow(Z_FIN_AST - curPos.z, 2);
        individuals[threadId].velDiff =  pow(VR_FIN_AST - curPos.vr, 2) + pow(VTHETA_FIN_AST - curPos.vtheta, 2) + pow(VZ_FIN_AST - curPos.vz, 2);


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

__host__ void initializePosition(Individual *individuals, int size){
    for(int i=0; i<size ;i++){
        individuals[i].initialize();
    }
}