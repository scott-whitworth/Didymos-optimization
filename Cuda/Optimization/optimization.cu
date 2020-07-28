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

// Used to see if the best individual is changing when compared to a previous individual. 
// Returns true if the currentBest is not equal to previousBest within a distinguishable difference
bool changeInBest(double previousBestPos, double previousBestVel, Individual currentBest, double distinguishRate) {
    //truncate is used here to compare floats via the distinguguishRate, to ensure that there has been relatively no change.
    if (trunc(previousBestPos/distinguishRate) != trunc(currentBest.posDiff/distinguishRate)) {
        return true;
    }
    else {
        /*
        if (trunc(previousBestVel/distinguishRate) != trunc(currentBest.velDiff/distinguishRate)) {
            return true;
        }
        else return false;
        */
        return false;
    }
}

// Assumes pool is sorted array of Individuals, used in determining if the loop continues
// Output: Returns true if top ten individuals within the pool are within the tolerance
bool allWithinTolerance(double tolerance, Individual * pool, unsigned int currentGeneration, const cudaConstants* cConstants) {
    // Uses for loop to pinpoint which individual is not in tolerance and display it to the terminal
    for (int i = 0; i < cConstants->best_count; i++) {
        if(pool[i].getPosDiff(cConstants) >= tolerance) {  // This isn't ideal, Change to getCost once getCost gets fleshed out //if (pool[i].getCost() >= tolerance ) {
            return false;
        }
    }
    // If iterated through and all were within tolerance, success
    return true;
    // std::cout << "\nallWithinTolerance() returned a best posDiff of " << pool[0].posDiff << std::endl;
}

// The function that starts up and runs the genetic algorithm with a continous loop until the critera is met (number of individuals equal to best_count is below the threshold value)
double optimize(const cudaConstants* cConstants) {
    thruster<double> thrust(cConstants);
    double calcPerS = 0;

    time_t timeSeed = cConstants->time_seed;
    std::mt19937_64 rng(timeSeed);
    std::cout << "----------------------------------------------------------------------------------------------------" << std::endl;
       
    // Initialize the recording files if in record mode
    if (cConstants->record_mode == true) {
        initializeRecord(cConstants);
    }
     
    // input parameters for rk4Simple which are the same for each thread
    double timeInitial = 0; // the starting time of the trip is always defined as zero   
    double absTol = cConstants->rk_tol; // the tolerance is a constant number that is shared amongst all runs
    double stepSize = (orbitalPeriod - timeInitial) / cConstants->max_numsteps; // the starting step size- same for each run- note that the current step size varies throughout each run

    double currentAnneal = cConstants->anneal_initial;

    Individual *inputParameters = new Individual[cConstants->num_individuals]; // contains all input parameters besides those which are always common amongst every thread

    double previousBestPos = 0; // set to zero to ensure there is a difference between previousBest and currentBest on generation zero (see changeInBest function)
    double previousBestVel = 0;

    if (cConstants->random_start) {
        // Sets inputParameters to hold parameters that are randomly generated within a reasonable range
        for (int i = 0; i < cConstants->num_individuals; i++) { 
            inputParameters[i] = Individual(randomParameters(rng, cConstants), cConstants);
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
        // arrayCPU needs to be updated to handle the fact that OPTIM_VARS may be flexible
        double arrayCPU[numStarts][OPTIM_VARS];
        
        for (int i = 0; i < OPTIM_VARS; i++) { // rows
            for (int j = 0; j < numStarts; j++) { // columns
                starts.read( reinterpret_cast<char*>( &startDoubles ), sizeof startDoubles );
                arrayCPU[j][i] = startDoubles;
            }
        }
        starts.close();

         // set every thread's input parameters to a set of final values from CPU calculations for use as a good starting point
        for (int i = 0; i < cConstants->num_individuals; i++) {
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

            inputParameters[i] = Individual(example, cConstants);
        }
    }


    Individual *survivors = new Individual[cConstants->survivor_count]; // stores the winners of the head-to-head competition
    int newInd = cConstants->num_individuals; // the whole population is new the first time through the loop

    double generation = 0;    // A counter for number of generations calculated
    
    double currentDistance; // Contains value for how far away the best individual is from the tolerance value
    double tolerance = cConstants->pos_threshold; // Tolerance for what is an acceptable solution (currently just the position threshold which is furthest distance from the target allowed)
                                                  // This could eventually take into account velocity too and become a more complex calculation
    //dRate is used to help check for a change in anneal
    double dRate = 1.0e-8;

    bool convergence = false;

    // A do-while loop that continues until it is determined that the pool of inputParameters has reached desired tolerance level for enough individuals (best_count)
    do {
        callRK(newInd, cConstants->thread_block_size, inputParameters + (cConstants->num_individuals - newInd), timeInitial, stepSize, absTol, calcPerS, thrust, cConstants); // calculate trajectories for new individuals

        // if we got bad results reset the Individual to random starting values (it may still be used for crossover) and set the final position to be way off so it gets replaced by a new Individual
        for (int k = 0; k < cConstants->num_individuals; k++) {
            if (isnan(inputParameters[k].finalPos.r) || isnan(inputParameters[k].finalPos.theta) || isnan(inputParameters[k].finalPos.z) || isnan(inputParameters[k].finalPos.vr) || isnan(inputParameters[k].finalPos.vtheta) || isnan(inputParameters[k].finalPos.vz)) {
                std::cout << std::endl << std::endl << "NAN FOUND" << std::endl << std::endl;
                inputParameters[k] = Individual(randomParameters(rng, cConstants), cConstants);
                // Set to be a bad individual
                inputParameters[k].posDiff = 1.0;
                inputParameters[k].velDiff = 0.0;
             }
            // calculate its new cost function
            inputParameters[k].getCost(cConstants);
        }
        // Note to future development, should shuffle and sort be within selectWinners method?
        std::shuffle(inputParameters, inputParameters + cConstants->num_individuals, rng); // shuffle the Individiuals to use random members for the competition
        selectSurvivors(inputParameters, cConstants->survivor_count, survivors); // Choose which individuals are in survivors, not necessarrily only the best ones
        std::sort(inputParameters, inputParameters + cConstants->num_individuals); // put the individuals in order so we can replace the worst ones

        // Display a '.' to the terminal to show that a generation has been performed
        // This also serves to visually seperate the generation display on the terminal screen
        std::cout << '.';
        // std::cout << "\nGeneration " << static_cast<int>(generation) << " returned a best posDiff of " << inputParameters[0].posDiff << std::endl;


        // Calculate how far the pool is from the ideal cost value (currently is the positionalDifference of the best individual)
        currentDistance = inputParameters[0].posDiff; // Change this later to take into account more than just the best individual and its position difference

        double new_anneal = currentAnneal * (1 - tolerance / currentDistance);
        
        Individual currentBest;
        if (static_cast<int>(generation) % cConstants->change_check == 0) { // Compare current best individual to that from CHANGE_CHECK many generations ago. If they are the same, change size of mutations
            currentBest = inputParameters[0];
            //checks for Anneal to change
            if ( !(changeInBest(previousBestPos, previousBestVel, currentBest, dRate)) ) { // previousBest starts at 0 to ensure changeInBest = true on generation 0
                //this ensures that changeInBest never compares two zeros, thus keeping dRate in relevance as the posDiff lowers.   
                if (trunc(currentBest.posDiff/dRate) == 0) { 
                    while (trunc(currentBest.posDiff/dRate) == 0) {
                        dRate = dRate/10; 
                    }
                    std::cout << "\nnew dRate: " << dRate << std::endl;
                }
                currentAnneal = currentAnneal * cConstants->anneal_factor;
                std::cout << "\nnew anneal: " << currentAnneal << std::endl;
            }
            previousBestPos = currentBest.posDiff;
            previousBestVel = currentBest.velDiff;
        }

        // If in recording mode and write_freq reached, call the record method
        if (static_cast<int>(generation) % cConstants->write_freq == 0 && cConstants->record_mode == true) {
            recordGenerationPerformance(cConstants, inputParameters, generation, new_anneal, cConstants->num_individuals, thrust);
        }
        // // Only call terminalDisplay every DISP_FREQ, not every single generation
        if ( static_cast<int>(generation) % cConstants->disp_freq == 0) {
            terminalDisplay(inputParameters[0], generation);
        }

        // Before replacing new individuals, determine whether all are within tolerance
        convergence = allWithinTolerance(tolerance, inputParameters, generation, cConstants);

        // Create a new generation and increment the generation counter
        newInd = newGeneration(survivors, inputParameters, cConstants->survivor_count, cConstants->num_individuals, new_anneal, cConstants, thrust, rng, generation);
        ++generation;
        
        // If the current distance is still higher than the tolerance we find acceptable, perform the loop again
    } while ( !convergence && generation < cConstants->max_generations);

    // Only call finalRecord if the results actually have reached the threshold
    if (allWithinTolerance(tolerance, inputParameters, generation, cConstants)) {
        finalRecord(cConstants, inputParameters, static_cast<int>(generation), thrust);
    }

    delete [] inputParameters;
    delete [] survivors;

    return calcPerS;
}

int main () {
    // display GPU properties and ensure we are using the right one
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    std::cout << "\n\nDevice Number: 0 \n";
    std::cout << "- Device name: " << prop.name << std::endl << std::endl;
    cudaSetDevice(0);
    
    cudaConstants * cConstants = new cudaConstants("../Config_Constants/genetic.config"); // Declare the genetic constants used, with file path being used to receive initial values

    double zero_seed = cConstants->time_seed;
    // Perform the optimization with optimize function
    for (int run = 0; run < cConstants->run_count; run++) {
            // Adjust the time_seed so it is unique based on each run
            cConstants->time_seed = zero_seed + run*100;

            // Display contents of cConstants being used for this run and how many runs
            std::cout << *cConstants << std::endl;
            std::cout << "Performing run #" << run << " with current parameters\n";

            // pre-calculate a table of Earth's position within possible mission time range
            launchCon = new EarthInfo(cConstants); // a global variable to hold Earth's position over time, assigned new info for this run of parameters
        
            // File stream for outputting values that were calculated in EarthInfo constructor
            /*if (cConstants->record_mode == true) {
                recordEarthData(cConstants, run);
            }*/
            optimize(cConstants);

            delete launchCon; // Deallocate launchCon info for this run as it may be needing a different time range in the next run
    }

    // Now that the optimize function is done (assumed that optimize() also records it), deallocate memory of the cudaConstants
    delete cConstants;
    
    return 0;
}