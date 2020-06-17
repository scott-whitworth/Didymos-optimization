// Didymos Optimization Project using CUDA and a genetic algorithm

#include "../constants.h" //used for wetMass
#include "../Earth_calculations/orbitalMotion.h" //used for trajectory() and trajectoryPrint()
#include "../Earth_calculations/earthInfo.h"
#include "../Runge_Kutta/runge_kuttaCUDA.cuh" //for testing rk4simple
#include "../Config_Constants/config.h"

#include <iostream> // cout
#include <iomanip> //used for setw(), sets spaces between values output
#include <time.h> //for seeding the random number generator
#include <random>
#include <chrono>

//#define SURVIVOR_COUNT 240 // number of individuals to use for crossover each generation--MUST BE DIVISIBLE BY 2 TO PAIR OFF FOR CROSSOVER
// 240 (survivors) / 2 (parents per pair) * 8 (offspring per pair) = 960 = half of 1920 --for k620 GPU
#define SURVIVOR_COUNT 360 // number of individuals to use for crossover each generation--MUST BE DIVISIBLE BY 2 TO PA

// Used to see if the best individual is changing
// Returns true if the currentBest is not equal to previousBest
bool changeInBest(double previousBest, double currentBest) {
    if (previousBest != currentBest) {
        return true;
    }
    else {
        return false;
    }
}

// Utility function to observe the trend of best individual in the algorithm through the generations
// Input: Two ofstreams (one to .csv file and another to binary), current generation number, best individual, and annealing value derived to be used in next generation crossover/mutation
// Output: The two streams are appended the individual's information and anneal value
void writeCurrentBestToFile(std::ofstream& ExcelOutput, std::ofstream& BinOutput, unsigned int &currentGeneration, Individual &individual, double& annealing ) {
    // Output the information to excel spreadsheet
    ExcelOutput << currentGeneration << ','
                << individual.posDiff << ',' << individual.velDiff << ',' // The positional and velocity difference
                << individual.finalPos.r << ',' << individual.finalPos.theta << ',' << individual.finalPos.z << ',' // Final position
                << individual.finalPos.vr << ',' << individual.finalPos.vtheta << ',' << individual.finalPos.vz << ',' // Final velocity
                << individual.startParams.y0.r << ',' << individual.startParams.y0.theta << ',' << individual.startParams.y0.z << ',' // Starting position
                << individual.startParams.y0.vr << ',' << individual.startParams.y0.vtheta << ',' << individual.startParams.y0.vz << ',' // Starting velocity
                << individual.startParams.alpha << ',' << individual.startParams.beta << ',' << individual.startParams.zeta << ',' // alpha, beta, zeta
                << annealing << "," << individual.startParams.tripTime << std::endl; // Annealing value for next generation and triptime (in that order to maintain continuity with bin file)
 
    // Output the information to binary file for use in the MATLAB code, line breaks and spaces added to help with readibility
    BinOutput.write( (char*)& currentGeneration, sizeof(unsigned int));
    // posDiff and velDiff
    BinOutput.write( (char*)& individual.posDiff, sizeof(double));
    BinOutput.write( (char*)& individual.velDiff, sizeof(double));
    // Position and velocity information
    BinOutput.write( (char*)& individual.finalPos.r,            sizeof(double));
    BinOutput.write( (char*)& individual.finalPos.theta,        sizeof(double));
    BinOutput.write( (char*)& individual.finalPos.z,            sizeof(double));
    BinOutput.write( (char*)& individual.finalPos.vr,           sizeof(double));
    BinOutput.write( (char*)& individual.finalPos.vtheta,       sizeof(double));
    BinOutput.write( (char*)& individual.finalPos.vz,           sizeof(double));
    BinOutput.write( (char*)& individual.startParams.y0.r,      sizeof(double));
    BinOutput.write( (char*)& individual.startParams.y0.theta,  sizeof(double));
    BinOutput.write( (char*)& individual.startParams.y0.z,      sizeof(double));
    BinOutput.write( (char*)& individual.startParams.y0.vr,     sizeof(double));
    BinOutput.write( (char*)& individual.startParams.y0.vtheta, sizeof(double));
    BinOutput.write( (char*)& individual.startParams.y0.vz,     sizeof(double));
    // Alpha, Beta, Zeta, Annealing, Triptime
    BinOutput.write( (char*)& individual.startParams.alpha,  sizeof(double));
    BinOutput.write( (char*)& individual.startParams.beta,   sizeof(double));
    BinOutput.write( (char*)& individual.startParams.zeta,   sizeof(double));
    BinOutput.write((char*)& annealing, sizeof(double));
    BinOutput.write((char*)& individual.startParams.tripTime, sizeof(double));
}

// Utility function to display the currently best individual onto the terminal while the algorithm is still running
// Input: Individual to be displayed (assumed to be the best individual of the pool) and the value for the current generation iterated
void terminalDisplay(Individual& individual, unsigned int currentGeneration) {
    std::cout << "\nGeneration: " << currentGeneration << std::endl;
    std::cout << "Best individual:" << std::endl;
    std::cout << "\tposDiff: " << individual.posDiff << std::endl;
    std::cout << "\tvelDiff: " << individual.velDiff << std::endl;
    std::cout << "\tcost: "    << individual.getCost() << std::endl;
}

// Assumes pool is sorted array of Individuals, used in determining if the loop continues
// Output: Returns true if top ten individuals within the pool are within the tolerance
bool allWithinTolerance(double tolerance, Individual * pool, unsigned int currentGeneration, geneticConstants& gConstant) {
    // Uses for loop to pinpoint which individual is not in tolerance and display it to the terminal
    for (int i = 0; i < gConstant.best_count; i++) {
        if (pool[i].getCost() >= tolerance ) {
            return false;
        }
    }
    // If iterated through and all were within tolerance, success
    return true;
}

// The function that starts up and runs the genetic algorithm with a continous loop until the critera is met (number of individuals equal to best_count is below the threshold value)
double optimize(const int numThreads, const int blockThreads, geneticConstants& gConstant) {
    double calcPerS = 0;

    time_t timeSeed = gConstant.time_seed;
    

    std::cout << "Seed for this run: " << timeSeed << std::endl; // note there are other mt_rands in the code that use different seeds
    std::cout << "------------------------------------------------------------------------" << std::endl;
    std::mt19937_64 mt_rand(timeSeed);

    // input parameters for rk4Simple which are the same for each thread
    double timeInitial = 0; // the starting time of the trip is always defined as zero   
    double absTol = RK_TOL; // the tolerance is a constant number that is shared amongst all runs
    double stepSize = (orbitalPeriod - timeInitial) / MAX_NUMSTEPS; // the starting step size- same for each run- note that the current step size varies throughout each run

    double currentAnneal = gConstant.anneal_initial;

    Individual *inputParameters = new Individual[numThreads]; // contains all input parameters besides those which are always common amongst every thread

    const int numStarts = 14; // the number of different sets of starting parameters in the input file

    std::ifstream starts;
    starts.open("../optimizedVector.bin", std::ifstream::in|std::ios::binary); // a file containing the final parameters of converged results from CPU calculations

    
    double previousBest = 0; // set to zero to ensure there is a difference between previousBest and currentBest on generation zero (see changeInBest function)

    // sort the data into 2 dimensions
    // one row is one set of starting parameters
    // each column is a specific variable:
    //    0-6 gamma
    //    7-9 tau
    //    10-12 launch angles
    //    13 trip time
    //    14-19 coast
    double startDoubles;
    double arrayCPU[numStarts][OPTIM_VARS];
    for (int i = 0; i < OPTIM_VARS; i++) { // rows
        for (int j = 0; j < numStarts; j++) { // columns
            starts.read( reinterpret_cast<char*>( &startDoubles ), sizeof startDoubles );
            arrayCPU[j][i] = startDoubles;
        }
    }
    starts.close();

     // set every thread's input parameters to a set of final values from CPU calculations for use as a good starting point
    for (int i = 0; i < numThreads; i++) {
        int row = mt_rand() % numStarts; // Choose a random row to get the parameters from

        double tripTime = arrayCPU[row][TRIPTIME_OFFSET];
        double alpha = arrayCPU[row][ALPHA_OFFSET];
        double beta = arrayCPU[row][BETA_OFFSET];
        double zeta = arrayCPU[row][ZETA_OFFSET];

        coefficients<double> testcoeff;
        for (int j = 0; j < testcoeff.gammaSize; j++) {
            testcoeff.gamma[j] = arrayCPU[row][j];
        }

        for (int j = 0; j < testcoeff.tauSize; j++) {
            testcoeff.tau[j] =  arrayCPU[row][j+7];
        }

        for (int j = 0; j < testcoeff.coastSize; j++) {
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

    // setup output of generation results over time onto a .csv file
    std::ofstream generationPerformanceExcel;
    generationPerformanceExcel.open("BestInGenerations.csv");
    // Set first row in the file be a header for the columns
    generationPerformanceExcel << "Gen #" << "," << "posDiff" << "," << "velDiff" << "," 
                               << "rFinal" << "," << "thetaFinal" << "," << "zFinal" << "," << "vrFinal" << "," << "vthetaFinal" << "," << "vzFinal" << ","
                               << "rInitial" << "," << "thetaInitial" << "," << "zInitial" << ","<< "vrInitial" << "," << "vthetaInitial" << "," << "vzInitial" << ","
                               << "alpha" << "," << "beta" << "," << "zeta" << "," << "anneal" << "," << "tripTime" << "\n";

    // setup output of generation results over time onto a .bin file
    std::ofstream generationPerformanceBin("BestInGenerations.bin", std::ios::binary);

    unsigned int generation = 0;    // A counter for number of generations calculated
    
    // A do-while loop that continues until it is determined that the pool of inputParameters has reached desired tolerance level for enough individuals (best_count)
    
    double currentDistance; // Contains value for how far away the best individual is from the tolerance value
    double tolerance = gConstant.pos_threshold; // Tolerance for what is an acceptable solution (currently just the position threshold which is furthest distance from the target allowed)
                                                // This could eventually take into account velocity too and become a more complex calculation

    do { // Set as a do while loop so that the algorithm is set to run atleast once
        // initialize positions for the new individuals starting at the index of the first new one and going to the end of the array
        initializePosition(inputParameters + (numThreads - newInd), newInd);

        callRK(newInd, blockThreads, inputParameters + (numThreads - newInd), timeInitial, stepSize, absTol, calcPerS); // calculate trajectories for new individuals

        // if we got bad results reset the Individual to random starting values (it may still be used for crossover) and set the final position to be way off so it gets replaced by a new Individual
        for (int k = 0; k < numThreads; k++) { 
            if (isnan(inputParameters[k].finalPos.r) || isnan(inputParameters[k].finalPos.theta) || isnan(inputParameters[k].finalPos.z) 
                 || isnan(inputParameters[k].finalPos.vr) || isnan(inputParameters[k].finalPos.vtheta) || isnan(inputParameters[k].finalPos.vz)){
                
                std::cout << std::endl << std::endl << "NAN FOUND" << std::endl << std::endl;

                double tripTime = 365*24*3600*(std::rand() % 10001 / 10000.0 + 1.0);
                double alpha = (mt_rand() % 629) / 100.0 - 3.14;
                double beta = (mt_rand() % 629) / 100.0 - 3.14;
                double zeta = (mt_rand() % 315) / 100.0 - 1.57;
        
                coefficients<double> testcoeff;
                for (int j = 0; j < testcoeff.gammaSize; j++) {
                    testcoeff.gamma[j] = mt_rand() % 201/10.0 - 10.0;
                }
                for (int j = 0; j < testcoeff.tauSize; j++) {
                    testcoeff.tau[j] = mt_rand() % 201/10.0 - 10.0;
                }
                for (int j = 0; j < testcoeff.coastSize; j++) {
                    testcoeff.coast[j] = mt_rand() % 201/10.0 - 10.0;
                }
            
                rkParameters<double> example(tripTime, alpha, beta, zeta, testcoeff); 
        
                inputParameters[k].startParams = example;

                inputParameters[k].posDiff = 1.0e10;
                inputParameters[k].velDiff = 0.0;
             }
        }
        // Note to future development, should shuffle and sort be within selectWinners method?
        std::shuffle(inputParameters, inputParameters + numThreads, mt_rand); // shuffle the Individiuals to use random members for the competition
        selectWinners(inputParameters, SURVIVOR_COUNT, survivors); // Choose which individuals are in survivors, not necessarrily only the best ones
        std::sort(inputParameters, inputParameters + numThreads); // put the individuals in order so we can replace the worst ones

        // Calculate how far the pool is from the ideal cost value (currently is the positionalDifference of the best individual)
        currentDistance = inputParameters[0].posDiff; // Change this later to take into account more than just the best individual and its position difference

        // the annealing rate passed in is scaled between ANNEAL_MAX and ANNEAL_MIN, dependent on the ratio between the tolerance and current distance from the tolerance
        // annealMax and annealMin change from the initial ANNEAL_MAX and ANNEAL_MIN whenever CHANGE_CHECK many generations pass without changing the best individual
        
        // double new_anneal =  annealMax - tolerance / currentDistance * (annealMax - annealMin); old way of calculating new_anneal that was made simpler
        double new_anneal = currentAnneal * (1 - tolerance / currentDistance);
        
        double currentBest;
        if (generation % gConstant.change_check == 0) { // Compare current best individual to that from CHANGE_CHECK many generations ago. If they are the same, change size of mutations
            currentBest = inputParameters[0].posDiff;
            if ( !(changeInBest(previousBest, currentBest)) ) { // previousBest starts at 0 to ensure changeInBest = true on generation 0
                currentAnneal = currentAnneal * gConstant.anneal_factor;
            }
            previousBest = inputParameters[0].posDiff;
        }

        // Display a '.' to the terminal to show that a generation has been performed
        // This also serves to visually seperate the generation display on the terminal screen
        std::cout << '.';

        // Write the best and worst Individuals in every 1,000 generations into the files to view progress over generations
        if (generation % gConstant.write_freq == 0) {
            writeCurrentBestToFile(generationPerformanceExcel, generationPerformanceBin, generation, inputParameters[0], new_anneal);
        }

        // Only call terminalDisplay every DISP_FREQ, not every single generation
        if ( generation % gConstant.disp_freq == 0) {
            terminalDisplay(inputParameters[0], generation);
        }

        // Create a new generation and increment the generation counter
        newInd = crossover(survivors, inputParameters, SURVIVOR_COUNT, numThreads, new_anneal, gConstant);
        ++generation;
        
        // If the current distance is still higher than the tolerance we find acceptable, perform the loop again
    } while ( !allWithinTolerance(tolerance, inputParameters, generation, gConstant) );


    
    // output the best Individuals of the final generation, using writeTrajectoryToFile()
    // Files outputted allows plotting of solutions in matlab
    double *start = new double[OPTIM_VARS];
    double cost = 0;

    // Output to excel
    double annealPlacement = 0; //setting anneal to be a placeholder value that has no real meaning as there will be no next generation for anneal to impact
    writeCurrentBestToFile(generationPerformanceExcel, generationPerformanceBin, generation, inputParameters[0], annealPlacement);

    // Write the best individuals with best_count in total outputted in seperate binary files
    for (int i = 0; i < gConstant.best_count; i++) {
        for (int j = 0; j < inputParameters[i].startParams.coeff.gammaSize; j++) {
            start[GAMMA_OFFSET + j] = inputParameters[i].startParams.coeff.gamma[j];
        }
        for (int j = 0; j < inputParameters[i].startParams.coeff.tauSize; j++) {
            start[TAU_OFFSET + j] = inputParameters[i].startParams.coeff.tau[j];
        }
        for (int j = 0; j < inputParameters[i].startParams.coeff.coastSize; j++) {
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

    // Close the performance files now that the algorithm is finished
    generationPerformanceExcel.close();
    generationPerformanceBin.close();

    delete [] inputParameters;
    delete [] survivors;
    delete start;

    return calcPerS;
}

int main () {
    // display GPU properties and ensure we are using the right one
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    std::cout << "Device Number: 0 \n";
    std::cout << "- Device name: " << prop.name << std::endl;
    cudaSetDevice(0);
    
    // pre-calculate a table of Earth's position within possible mission time range
    //----------------------------------------------------------------
    // Define variables to be passed into EarthInfo
    double startTime = 15778800; // 0.5 year (s)
    double endTime = 78894000; // 2.5 years (s)
    double timeRes = 3600; // (s) position of earth is calculated for every hour

    launchCon = new EarthInfo(startTime, endTime, timeRes); // a global variable to hold Earth's position over time
    //----------------------------------------------------------------
    // Define the number of threads/individuals that will be used in optimize
    int blockThreads = 32;
    int numThreads = 2880; // the number of cores on a Tesla k40
    //int numThreads = 1920; // 384 cores on K620 * 5 = 1920

    // Initialize the type of thruster to be used
    thruster<double> thrust(0);

    //std::ofstream efficiencyGraph; // for viewing how many runge-kuttas ran per second for each combination of threads per block and total threads 
    //efficiencyGraph.open("efficiencyGraph.csv");
    std::cout << std::endl << "running optimize() with " << blockThreads << " threads per block and " << numThreads << " total threads" << std::endl;
    
    geneticConstants gConstant("genetic.config"); // Declare the genetic constants used, with file path being used

    optimize(numThreads, blockThreads, gConstant);

    //efficiencyGraph << blockThreads << "," << numThreads << "," << calcPerS  << "\n";
    //efficiencyGraph.close();
    
    delete launchCon;

    return 0;
}

double optimize(const int numThreads, const int blockThreads, thruster<double> thrust) {
    double calcPerS = 0;
    time_t timeSeed = time(0); // 1234567890;
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
    for (int i = 0; i < OPTIM_VARS; i++) { // rows
        for (int j = 0; j < numStarts; j++) { // columns
            starts.read( reinterpret_cast<char*>( &startDoubles ), sizeof startDoubles );
            arrayCPU[j][i] = startDoubles;
        }
    }
    starts.close();

     // set every thread's input parameters to a set of final values from CPU calculations for use as a good starting point
    for (int i = 0; i < numThreads; i++) {
        int row = mt_rand() % numStarts;

        double tripTime = arrayCPU[row][13];

        double alpha = arrayCPU[row][10];

        double beta = arrayCPU[row][11];

        double zeta = arrayCPU[row][12];

        coefficients<double> testcoeff;

        if (thrust.type) {
            for (int j = 0; j < testcoeff.gammaSize; j++) {
                testcoeff.gamma[j] = arrayCPU[row][j];
            }

            for (int j = 0; j < testcoeff.tauSize; j++) {
                testcoeff.tau[j] =  arrayCPU[row][j+7];
            }

            for (int j = 0; j < testcoeff.coastSize; j++) {
                testcoeff.coast[j] = arrayCPU[row][j+14];
            }
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
    
    double posDiffRange = 0, velDiffRange = 0, prevBestPos = 0, prevBestVel = 0, prevWorstPos = 0, prevWorstVel = 0;

    //local variables for anneal
    double annealMax = ANNEAL_MAX;
    double annealMin = ANNEAL_MIN;
    double distinguishRate = 1.0e-7;

    // Initialize a generation counter and convergence flag
    double i = 1;
    bool convgFlag = 0;

    std::ofstream individualDifference;
    individualDifference.open("individualDifference.bin", std::ios::binary);

    while (!convgFlag) {
        // initialize positions for the new individuals starting at the index of the first new one and going to the end of the array
        initializePosition(inputParameters + (numThreads - newInd), newInd);

        callRK(newInd, blockThreads, inputParameters + (numThreads - newInd), timeInitial, stepSize, absTol, calcPerS, thrust); // calculate trajectories for new individuals

        for (int k = 0; k < numThreads; k++) { // if we got bad results reset the Individual to random starting values (it may still be used for crossover) 
                                            // and set the final position to be way off so it gets replaced by a new Individual
            if (isnan(inputParameters[k].finalPos.r) || isnan(inputParameters[k].finalPos.theta) || isnan(inputParameters[k].finalPos.z) 
                 || isnan(inputParameters[k].finalPos.vr) || isnan(inputParameters[k].finalPos.vtheta) || isnan(inputParameters[k].finalPos.vz)){
                
                std::cout << std::endl << std::endl << "NAN FOUND" << std::endl << std::endl;

                double tripTime = 365*24*3600*(std::rand() % 10001 / 10000.0 + 1.0); // (seconds In A Year) * (random decimal)
                double alpha = (mt_rand() % 629) / 100.0 - M_PI;
                double beta = (mt_rand() % 629) / 100.0 - M_PI;
                double zeta = (mt_rand() % 315) / 100.0 - M_PI / 2;
        
                coefficients<double> testcoeff;
                if (thrust.type) {
                    for (int j = 0; j < testcoeff.gammaSize; j++) {
                        testcoeff.gamma[j] = mt_rand() % 201/10.0 - 10.0;
                    }
                    for (int j = 0; j < testcoeff.tauSize; j++) {
                        testcoeff.tau[j] = mt_rand() % 201/10.0 - 10.0;
                    }
                    for (int j = 0; j < testcoeff.coastSize; j++) {
                        testcoeff.coast[j] = mt_rand() % 201/10.0 - 10.0;
                    }
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

        // Calculate the current generation's cost function range
        posDiffRange = posCost(inputParameters, numThreads);
        velDiffRange = velCost(inputParameters, numThreads);

        // Check for convergence
        convgFlag = posConverge(inputParameters);

        // finding the best variable to change in the best Individual
        // bestChange() TO BE USED HERE


        // Display a '.' to the terminal to show that a generation has been calculated
        // if it is not the 100th generation this serves to show that a generation was calculated and survivors selected
        // This also serves to visually seperate the generation display on the terminal screen
        std::cout << '.';

        double previousAnneal =  annealMax - static_cast<double>(i) / (generationsNum - 1) * (annealMax - annealMin);
        // Display and print Individuals' pos and vel difference every 100 generations to terminal and .csv file
        if (static_cast<int>(i) % 20 == 0) { 
            // Display the cost function range within every 100th generation
            
            std::cout << '\n';
            std::cout << "generation: " << i << std::endl;
            std::cout << "posDiffRange: " << posDiffRange << std::endl;
            std::cout << "velDiffRange: " << velDiffRange << std::endl;
            
            std::cout << "posDiffRange change over 20 gens: " << posDiffRange - abs(prevBestPos - prevWorstPos) <<std::endl;
            std::cout << "velDiffRange change over 20 gens: " << velDiffRange - abs(prevBestVel - prevWorstVel) <<std::endl;

            std::cout << "best:" << std::endl;
            std::cout << "\tposDiff: " << inputParameters[0].posDiff << std::endl;
            std::cout << "\tvelDiff: " << inputParameters[0].velDiff << std::endl;
            std::cout << "worst:" << std::endl;
            std::cout << "\tposDiff: " << inputParameters[numThreads - 1].posDiff << std::endl;
            std::cout << "\tvelDiff: " << inputParameters[numThreads - 1].velDiff << std::endl;

            
            if(static_cast<int>(i) % 50 == 0 && distinguishableDifference(prevBestPos, inputParameters[0].posDiff, distinguishRate)) {
                //half anneal  max and min
                annealMax = annealMax / ANNEAL_FACTOR;
                annealMin = annealMin / ANNEAL_FACTOR;
                if(trunc(inputParameters[0].posDiff/distinguishRate)==0) {
                    distinguishRate = distinguishRate/10;
                }
                std::cout << "New anneal max: " << annealMax << "  New anneal min: " << annealMin << " New distinguishRate: " << distinguishRate << std::endl;
            } 

            prevBestPos = inputParameters[0].posDiff;
            prevBestVel = inputParameters[0].velDiff;
            prevWorstPos = inputParameters[numThreads-1].posDiff;
            prevWorstVel = inputParameters[numThreads-1].velDiff;
        }
            // Append the best and worst Individuals into a bin file to view progress over generations
            writeProgressToFile(individualDifference, inputParameters, i, 0, previousAnneal);
            writeProgressToFile(individualDifference, inputParameters, i, numThreads-1, previousAnneal);

            // individualDifferenceCSV << i << "," << inputParameters[0].posDiff << "," << inputParameters[0].velDiff << "," << 
            //     inputParameters[0].finalPos.r << "," << inputParameters[0].finalPos.theta << "," << inputParameters[0].finalPos.z << "," << 
            //     inputParameters[0].finalPos.vr << "," << inputParameters[0].finalPos.vtheta << "," << inputParameters[0].finalPos.vz << "," << 
            //     inputParameters[0].startParams.y0.r << "," << inputParameters[0].startParams.y0.theta << "," << inputParameters[0].startParams.y0.z << "," << 
            //     inputParameters[0].startParams.y0.vr << "," << inputParameters[0].startParams.y0.vtheta << "," << inputParameters[0].startParams.y0.vz << "," << 
            //     inputParameters[0].startParams.alpha << "," << inputParameters[0].startParams.beta << "," << inputParameters[0].startParams.zeta << "," << 
            //     previousAnneal << "," << inputParameters[0].startParams.tripTime << "\n";

        // the annnealing rate passed in is scaled between ANNEAL_MAX and ANNEAL_MIN depending on which generation this is
        double new_anneal =  annealMax - static_cast<double>(i) / (generationsNum - 1) * (annealMax - annealMin);
       

        newInd = crossover(survivors, inputParameters, SURVIVOR_COUNT, numThreads, new_anneal, thrust);

        // Step into the next generation
        i++;
        
    }


    // output the best Individual of the final generation, using writeTrajectoryToFile()
    // Files outputted allows plotting of solutions in matlab
    double *start = new double[OPTIM_VARS];
    double cost = 0;

    for (int i = 0; i < 1; i++) { 
        if (thrust.type) {
            for (int j = 0; j < inputParameters[i].startParams.coeff.gammaSize; j++) {
                start[GAMMA_OFFSET + j] = inputParameters[i].startParams.coeff.gamma[j];
            }
            for (int j = 0; j < inputParameters[i].startParams.coeff.tauSize; j++) {
                start[TAU_OFFSET + j] = inputParameters[i].startParams.coeff.tau[j];
            }
            for (int j = 0; j < inputParameters[i].startParams.coeff.coastSize; j++) {
                start[COAST_OFFSET + j] = inputParameters[i].startParams.coeff.coast[j];
            }
        }

        start[TRIPTIME_OFFSET] = inputParameters[i].startParams.tripTime;
        start[ALPHA_OFFSET] = inputParameters[i].startParams.alpha;
        start[BETA_OFFSET] = inputParameters[i].startParams.beta;
        start[ZETA_OFFSET] = inputParameters[i].startParams.zeta;

        cost = inputParameters[i].posDiff; // just look at position difference here for now
        // could instead use a ratio between position and velocity differnce as done in comparison of Individuals
        writeTrajectoryToFile(start, cost, i + 1, thrust);
    }

    individualDifference.close();

    delete [] inputParameters;
    delete [] survivors;

    return calcPerS;
}
