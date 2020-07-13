// Didymos Optimization Project using CUDA and a genetic algorithm

#include "../Earth_calculations/earthInfo.h"
#include "../Runge_Kutta/runge_kuttaCUDA.cuh" //for testing rk4simple
#include "../Config_Constants/config.h"
#include "../output.h"

#include <iostream> // cout
#include <iomanip>  // used for setw(), sets spaces between values output
#include <time.h>   // for seeding the random number generator
#include <random>
#include <chrono>

#define SECONDS_IN_YEAR 365.25*24*3600

// Used to see if the best individual is changing
// Returns true if the currentBest is not equal to previousBest
bool changeInBest(double previousBestPos, double previousBestVel, Individual currentBest, double dRate) {
    if (trunc(previousBestPos/dRate) != trunc(currentBest.posDiff/dRate)) {
        return true;
    }
    else {
        if (trunc(previousBestVel/dRate) != trunc(currentBest.velDiff/dRate)) {
            return true;
        }
        else return false;
    }
}

// Utility function to display the currently best individual onto the terminal while the algorithm is still running
// Input: Individual to be displayed (assumed to be the best individual of the pool) and the value for the current generation iterated
void terminalDisplay(Individual& individual, unsigned int currentGeneration) {
    std::cout << "\nGeneration: " << currentGeneration << std::endl;
    std::cout << "Best individual:" << std::endl;
    std::cout << "\tposDiff: " << individual.posDiff << std::endl;
    std::cout << "\tvelDiff: " << individual.velDiff << std::endl;
    std::cout << "\tcost: "    << individual.cost << std::endl;
}

// Assumes pool is sorted array of Individuals, used in determining if the loop continues
// Output: Returns true if top ten individuals within the pool are within the tolerance
bool allWithinTolerance(double tolerance, Individual * pool, unsigned int currentGeneration, const cudaConstants* cConstants) {
    // Uses for loop to pinpoint which individual is not in tolerance and display it to the terminal
    for (int i = 0; i < cConstants->best_count; i++) {
        if(pool[i].posDiff >= cConstants->pos_threshold) {  // This isn't ideal, Change to getCost once getCost gets fleshed out //if (pool[i].getCost() >= tolerance ) {
            return false;
        }
    }
    // If iterated through and all were within tolerance, success
    return true;
}

// The function that starts up and runs the genetic algorithm with a continous loop until the critera is met (number of individuals equal to best_count is below the threshold value)
double optimize(const int numThreads, const int blockThreads, const cudaConstants* cConstants, thruster<double> thrust) {
    double calcPerS = 0;

    time_t timeSeed = cConstants->time_seed;
    std::mt19937_64 rng(timeSeed);
    std::cout << "----------------------------------------------------------------------------------------------------" << std::endl;
    
    //sets up mutate file
    setMutateFile(cConstants);
    
    // input parameters for rk4Simple which are the same for each thread
    double timeInitial = 0; // the starting time of the trip is always defined as zero   
    double absTol = cConstants->rk_tol; // the tolerance is a constant number that is shared amongst all runs
    double stepSize = (orbitalPeriod - timeInitial) / cConstants->max_numsteps; // the starting step size- same for each run- note that the current step size varies throughout each run

    double currentAnneal = cConstants->anneal_initial;

    Individual *inputParameters = new Individual[numThreads]; // contains all input parameters besides those which are always common amongst every thread

    double previousBestPos = 0; // set to zero to ensure there is a difference between previousBest and currentBest on generation zero (see changeInBest function)
    double previousBestVel = 0;

    if (cConstants->random_start) {
        // Sets inputParameters to hold parameters that are randomly generated within a reasonable range
        for (int i = 0; i < numThreads; i++) { 
            double tripTime = SECONDS_IN_YEAR*(rng() % 10001 / 10000.0 + 1.0); // (1 <-> 2 years) * SECONDS_IN_YEAR
            double alpha = M_PI * 2*((static_cast<double>(rng()) / rng.max()) - 0.5); // -PI <-> PI
            double beta  = M_PI * ((static_cast<double>(rng()) / rng.max())); // 0 <-> PI
            double zeta  = M_PI * ((static_cast<double>(rng()) / rng.max()) - 0.5); // -PI/2 <-> PI/2

            coefficients<double> testcoeff;
            for (int j = 0; j < testcoeff.gammaSize; j++) {
                testcoeff.gamma[j] = rng() % 201/10.0 - 10.0; // -10.0 <-> 10.0
            }
            for (int j = 0; j < testcoeff.tauSize; j++) {
                testcoeff.tau[j] = rng() % 201/10.0 - 10.0; // -10.0 <-> 10.0
            }
            for (int j = 0; j < testcoeff.coastSize; j++) {
                testcoeff.coast[j] = rng() % 201/10.0 - 10.0; // -10.0 <-> 10.0
            }
        
            rkParameters<double> example(tripTime, alpha, beta, zeta, testcoeff); 
        
            inputParameters[i].startParams = example;
        }
    }
    // If not a random start, read from file using cConstants initial_start_file_address to get path
    else {
        // Sets inputParameters to hold initial individuals based from file optimizedVector.bin

        const int numStarts = 14; // the number of different sets of starting parameters in the input file

        std::ifstream starts;
        starts.open(cConstants->initial_start_file_address, std::ifstream::in|std::ios::binary); // a file containing the final parameters of converged results from CPU calculations        

        // sort the data into 2 dimensions
        // one row is one set of starting parameters
        // each column is a specific variable:
        double startDoubles;
        // arrayCPU needs to be updated to handle the fact that OPTIM_VARS now is defined from cConstants
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
            int row = rng() % numStarts; // Choose a random row to get the parameters from

            double tripTime = arrayCPU[row][TRIPTIME_OFFSET];
            double alpha = arrayCPU[row][ALPHA_OFFSET];
            double beta = arrayCPU[row][BETA_OFFSET];
            double zeta = arrayCPU[row][ZETA_OFFSET];

            coefficients<double> testcoeff;
            for (int j = 0; j < testcoeff.gammaSize; j++) {
                testcoeff.gamma[j] = arrayCPU[row][j + GAMMA_OFFSET];
            }

            for (int j = 0; j < testcoeff.tauSize; j++) {
                testcoeff.tau[j] =  arrayCPU[row][j + TAU_OFFSET];
            }

            for (int j = 0; j < testcoeff.coastSize; j++) {
                testcoeff.coast[j] = arrayCPU[row][j + COAST_OFFSET];
            }

            rkParameters<double> example(tripTime, alpha, beta, zeta, testcoeff); 

            inputParameters[i].startParams = example;
        }
    }


    Individual *survivors = new Individual[cConstants->survivor_count]; // stores the winners of the head-to-head competition
    int newInd = numThreads; // the whole population is new the first time through the loop

    // Initialize the recording files if in record mode
    if (cConstants->record_mode == true) {
        initializeRecord(cConstants);
    }

    double generation = 0;    // A counter for number of generations calculated
    
    // A do-while loop that continues until it is determined that the pool of inputParameters has reached desired tolerance level for enough individuals (best_count)
    
    double currentDistance; // Contains value for how far away the best individual is from the tolerance value
    double tolerance = cConstants->pos_threshold; // Tolerance for what is an acceptable solution (currently just the position threshold which is furthest distance from the target allowed)
                                                  // This could eventually take into account velocity too and become a more complex calculation
    double dRate = 1.0e-8;

    do { // Set as a do while loop so that the algorithm is set to run atleast once
        // initialize positions for the new individuals starting at the index of the first new one and going to the end of the array
        initializePosition(inputParameters + (numThreads - newInd), newInd, cConstants);
        callRK(newInd, blockThreads, inputParameters + (numThreads - newInd), timeInitial, stepSize, absTol, calcPerS, thrust, cConstants); // calculate trajectories for new individuals

        // if we got bad results reset the Individual to random starting values (it may still be used for crossover) and set the final position to be way off so it gets replaced by a new Individual
        for (int k = 0; k < numThreads; k++) {
            if (isnan(inputParameters[k].finalPos.r) || isnan(inputParameters[k].finalPos.theta) || isnan(inputParameters[k].finalPos.z) 
                 || isnan(inputParameters[k].finalPos.vr) || isnan(inputParameters[k].finalPos.vtheta) || isnan(inputParameters[k].finalPos.vz)){
                
                std::cout << std::endl << std::endl << "NAN FOUND" << std::endl << std::endl;

                double tripTime = SECONDS_IN_YEAR*(std::rand() % 10001 / 10000.0 + 1.0);
                double alpha = M_PI * 2*((static_cast<double>(rng()) / rng.max()) - 0.5); // -PI <-> PI
                double beta  = M_PI *   ((static_cast<double>(rng()) / rng.max())); // 0 <-> PI
                double zeta  = M_PI *   ((static_cast<double>(rng()) / rng.max()) - 0.5); // -PI/2 <-> PI/2

                coefficients<double> testcoeff;
                if (thrust.type) {
                    for (int j = 0; j < testcoeff.gammaSize; j++) {
                        testcoeff.gamma[j] = rng() % 201/10.0 - 10.0;
                    }
                    for (int j = 0; j < testcoeff.tauSize; j++) {
                        testcoeff.tau[j] = rng() % 201/10.0 - 10.0;
                    }
                    for (int j = 0; j < testcoeff.coastSize; j++) {
                        testcoeff.coast[j] = rng() % 201/10.0 - 10.0;
                    }
                }
            
                rkParameters<double> example(tripTime, alpha, beta, zeta, testcoeff); 
        
                inputParameters[k].startParams = example;

                // Set to be a bad individual
                inputParameters[k].posDiff = 1.0;
                inputParameters[k].velDiff = 0.0;
             }

            // calculate its new cost function
            inputParameters[k].getCost(cConstants);
        }

        // Note to future development, should shuffle and sort be within selectWinners method?
        std::shuffle(inputParameters, inputParameters + numThreads, rng); // shuffle the Individiuals to use random members for the competition
        selectSurvivors(inputParameters, cConstants->survivor_count, survivors); // Choose which individuals are in survivors, not necessarrily only the best ones
        std::sort(inputParameters, inputParameters + numThreads); // put the individuals in order so we can replace the worst ones

        // Display a '.' to the terminal to show that a generation has been performed
        // This also serves to visually seperate the generation display on the terminal screen
        std::cout << '.';

        // Calculate how far the pool is from the ideal cost value (currently is the positionalDifference of the best individual)
        currentDistance = inputParameters[0].posDiff; // Change this later to take into account more than just the best individual and its position difference

        double new_anneal = currentAnneal * (1 - tolerance / currentDistance);
        
        Individual currentBest;
        if (static_cast<int>(generation) % cConstants->change_check == 0) { // Compare current best individual to that from CHANGE_CHECK many generations ago. If they are the same, change size of mutations
            currentBest = inputParameters[0];
          
            if ( !(changeInBest(previousBestPos, previousBestVel, currentBest, dRate)) ) { // previousBest starts at 0 to ensure changeInBest = true on generation 0
                currentAnneal = currentAnneal * cConstants->anneal_factor;
                std::cout << "\nnew anneal: " << currentAnneal << std::endl;
                if(trunc(inputParameters[0].posDiff/dRate)==0) { dRate = dRate/10; }
            }
            previousBestPos = currentBest.posDiff;
            previousBestVel = currentBest.velDiff;
        }

        // If in recording mode and write_freq reached, call the record method
        if (static_cast<int>(generation) % cConstants->write_freq == 0 && cConstants->record_mode == true) {
            recordGenerationPerformance(cConstants, inputParameters, generation, new_anneal, numThreads);
        }

        // Only call terminalDisplay every DISP_FREQ, not every single generation
        if ( static_cast<int>(generation) % cConstants->disp_freq == 0) {
            terminalDisplay(inputParameters[0], generation);
        }

        // Create a new generation and increment the generation counter
        newInd = newGeneration(survivors, inputParameters, cConstants->survivor_count, numThreads, new_anneal, cConstants, thrust, rng, generation);
        ++generation;
        
        // If the current distance is still higher than the tolerance we find acceptable, perform the loop again
    } while ( !allWithinTolerance(tolerance, inputParameters, generation, cConstants) );

    
    // output the best Individuals of the final generation, using writeTrajectoryToFile()
    // Files outputted allows plotting of solutions in matlab
    double *start = new double[OPTIM_VARS];

    // Write the final best and worst performing individuals to their respective files    
    recordGenerationPerformance(cConstants, inputParameters, generation, 0, numThreads);
    
    // std::ofstream progressiveOutput;
    // progressiveOutput.open("progressiveAnalysis.csv", std::ios::app);
    // progressiveOutput << std::endl << "seed:," << cConstants->time_seed << ",  ,generations:," << static_cast<int>(generation) << std::endl;
    // progressiveOutput << "rank,posDiff (au),velDiff (au/s),tripTime (s),alpha (rad),beta (rad),zeta (rad),";
    // if (thrust.type) {
    //     progressiveOutput << "gamma_a0,gamma_a1,gamma_b1,gamme_a2,gamme_b2,gamma_a3,gamma_b3,";
    //     progressiveOutput << "tau_a0,tau_a1,tau_b1,";
    //     progressiveOutput << "coast_a0,coast_a1,coast_b1,coast_a2,coast_b2,";
    // }
    // progressiveOutput << std::endl;
    // // Write the best individuals with best_count in total outputted in seperate binary files
    // for (int i = 0; i < cConstants->best_count; i++) {
    //     for (int j = 0; j < inputParameters[i].startParams.coeff.gammaSize; j++) {
    //         start[GAMMA_OFFSET + j] = inputParameters[i].startParams.coeff.gamma[j];
    //     }
    //     for (int j = 0; j < inputParameters[i].startParams.coeff.tauSize; j++) {
    //         start[TAU_OFFSET + j] = inputParameters[i].startParams.coeff.tau[j];
    //     }
    //     for (int j = 0; j < inputParameters[i].startParams.coeff.coastSize; j++) {
    //         start[COAST_OFFSET + j] = inputParameters[i].startParams.coeff.coast[j];
    //     }

    //     start[TRIPTIME_OFFSET] = inputParameters[i].startParams.tripTime;
    //     start[ALPHA_OFFSET] = inputParameters[i].startParams.alpha;
    //     start[BETA_OFFSET] = inputParameters[i].startParams.beta;
    //     start[ZETA_OFFSET] = inputParameters[i].startParams.zeta;

    //     // could instead use a ratio between position and velocity differnce as done in comparison of Individuals
    //     writeTrajectoryToFile(start, i+1, thrust, cConstants);
    //     progressiveAnalysis(progressiveOutput,i+1,inputParameters[i],cConstants);
    // }
    // progressiveOutput << std::endl;
    // progressiveOutput.close();

    // Only output the final best individual
    for (int j = 0; j < inputParameters[0].startParams.coeff.gammaSize; j++) {
        start[GAMMA_OFFSET + j] = inputParameters[0].startParams.coeff.gamma[j];
    }
    for (int j = 0; j < inputParameters[0].startParams.coeff.tauSize; j++) {
        start[TAU_OFFSET + j] = inputParameters[0].startParams.coeff.tau[j];
    }
    for (int j = 0; j < inputParameters[0].startParams.coeff.coastSize; j++) {
        start[COAST_OFFSET + j] = inputParameters[0].startParams.coeff.coast[j];
    }

    start[TRIPTIME_OFFSET] = inputParameters[0].startParams.tripTime;
    start[ALPHA_OFFSET] = inputParameters[0].startParams.alpha;
    start[BETA_OFFSET] = inputParameters[0].startParams.beta;
    start[ZETA_OFFSET] = inputParameters[0].startParams.zeta;

    // could instead use a ratio between position and velocity differnce as done in comparison of Individuals
    writeTrajectoryToFile(start, 1, thrust, cConstants);

    delete [] inputParameters;
    delete [] survivors;
    delete start;

    return calcPerS;
}

int main () {
    // display GPU properties and ensure we are using the right one
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    std::cout << "\n\nDevice Number: 0 \n";
    std::cout << "- Device name: " << prop.name << std::endl << std::endl;
    cudaSetDevice(0);
    
    cudaConstants const * cConstants = new cudaConstants("../Config_Constants/genetic.config"); // Declare the genetic constants used, with file path being used
    // Display contents of cConstants resulting from reading the file onto the terminal
    std::cout << *cConstants << std::endl;

    // pre-calculate a table of Earth's position within possible mission time range
    //----------------------------------------------------------------
    // Define variables to be passed into EarthInfo that determines the range of time to be calculated, accessed from cConstants
    double startTime = cConstants->startTime;
    double endTime = cConstants->endTime; 
    double timeRes = cConstants->timeRes;

    launchCon = new EarthInfo(startTime, endTime, timeRes, cConstants); // a global variable to hold Earth's position over time

    // double timeStamp = startTime;
    
    // // File stream for outputting values that were calculated in EarthInfo constructor
    if (cConstants->record_mode) {
        std::ofstream earthValues;
        earthValues.open("EarthCheckValues.csv");
        // Set header row for the table to record values, with timeStamp
        earthValues << "TimeStamp, Radius, Theta, Z, vRadius, vTheta, vZ\n";
        while (timeStamp < endTime) {
            earthValues << timeStamp << "," << launchCon->getCondition(timeStamp);
            timeStamp += timeRes*24; // Increment to next day as timeRes is set to every hour
        }
        // Done recording earth calculations, close file and move on
        earthValues.close();
    }
    
    //----------------------------------------------------------------
    // Define the number of threads/individuals that will be used in optimize
    int blockThreads = cConstants->thread_block_size;
    int numThreads = cConstants->num_individuals;

    std::cout << std::endl << "running optimize() with " << blockThreads << " threads per block and " << numThreads << " total threads" << std::endl;

    thruster<double> thrust(cConstants);

    // Perform the optimization with optimize function
    optimize(numThreads, blockThreads, cConstants, thrust);

    // Now that the optimize function is done (assumed taht optimize() also records it), deallocate memory of the earth calculations and cudaConstants
    delete launchCon;
    delete cConstants;
    
    return 0;
}