// Didymos Optimization Project using CUDA and a genetic algorithm

#include "../Earth_calculations/orbitalMotion.h" //used for trajectory() and trajectoryPrint()
#include "../Earth_calculations/earthInfo.h"
#include "../Runge_Kutta/runge_kuttaCUDA.cuh" //for testing rk4simple
#include "../Config_Constants/config.h"

#include <iostream> // cout
#include <iomanip> //used for setw(), sets spaces between values output
#include <time.h> //for seeding the random number generator
#include <random>
#include <chrono>

#define PI 3.141592653

//#define SURVIVOR_COUNT 240 // number of individuals to use for crossover each generation--MUST BE DIVISIBLE BY 2 TO PAIR OFF FOR CROSSOVER
// 240 (survivors) / 2 (parents per pair) * 8 (offspring per pair) = 960 = half of 1920 --for k620 GPU
#define SURVIVOR_COUNT 360 // number of individuals to use for crossover each generation--MUST BE DIVISIBLE BY 2 TO PA
#define SECONDS_IN_YEAR 365*24*3600

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

// Utility function to observe the trend of best individual in the algorithm through the generations
// Input: Two ofstreams (one to .csv file and another to binary), current generation number, best individual, and annealing value derived to be used in next generation crossover/mutation
// Output: The two streams are appended the individual's information and anneal value
void writeIndividualToFiles(std::ofstream& ExcelOutput, std::ofstream& BinOutput, double &currentGeneration, Individual &individual, double& annealing ) {
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
    BinOutput.write( (char*)& currentGeneration, sizeof(double));
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

void writeThrustToFiles(std::ofstream& ExcelOutput, std::ofstream& BinOutput, double &currentGeneration, Individual &individual, const cudaConstants * cConstants) {
    ExcelOutput << currentGeneration << ',';
    for (int i = 0; i < GAMMA_ARRAY_SIZE; i++) {
        ExcelOutput << individual.startParams.coeff.gamma[i] << ',';
    }
    for (int i = 0; i < TAU_ARRAY_SIZE; i++) {
        ExcelOutput << individual.startParams.coeff.tau[i] << ',';
    }
    for (int i = 0; i < COAST_ARRAY_SIZE; i++) {
        ExcelOutput << individual.startParams.coeff.coast[i] << ',';
    }
    ExcelOutput << std::endl;

    BinOutput.write((char*)&currentGeneration, sizeof(double));
    for (int i = 0; i < GAMMA_ARRAY_SIZE; i++) {
        BinOutput.write((char*)&individual.startParams.coeff.gamma[i], sizeof(double));
    }
    for (int i = 0; i < TAU_ARRAY_SIZE; i++) {
        BinOutput.write((char*)&individual.startParams.coeff.tau[i], sizeof(double));
    }
    for (int i = 0; i < COAST_ARRAY_SIZE; i++) {
        BinOutput.write((char*)&individual.startParams.coeff.coast[i], sizeof(double));
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
        if(pool[i].posDiff >= cConstants->pos_threshold || pool[i].velDiff < cConstants->v_impact) {  // This isn't ideal, Change to getCost once getCost gets fleshed out //if (pool[i].getCost() >= tolerance ) {
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

    std::cout << "Seed for this run: " << timeSeed << std::endl; // note there are other mt_rands in the code that use different seeds
    std::cout << "------------------------------------------------------------------------" << std::endl;
    std::mt19937_64 mt_rand(timeSeed);

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
            double tripTime = SECONDS_IN_YEAR*(mt_rand() % 10001 / 10000.0 + 1.0); // (1 <-> 2 years) * SECONDS_IN_YEAR
            double alpha = PI * 2*((static_cast<double>(mt_rand()) / mt_rand.max()) - 0.5); // -PI <-> PI
            double beta  = PI * ((static_cast<double>(mt_rand()) / mt_rand.max())); // 0 <-> PI
            double zeta  = PI * ((static_cast<double>(mt_rand()) / mt_rand.max()) - 0.5); // -PI/2 <-> PI/2

            coefficients<double> testcoeff;
            for (int j = 0; j < testcoeff.gammaSize; j++) {
                testcoeff.gamma[j] = mt_rand() % 201/10.0 - 10.0; // -10.0 <-> 10.0
            }
            for (int j = 0; j < testcoeff.tauSize; j++) {
                testcoeff.tau[j] = mt_rand() % 201/10.0 - 10.0; // -10.0 <-> 10.0
            }
            for (int j = 0; j < testcoeff.coastSize; j++) {
                testcoeff.coast[j] = mt_rand() % 201/10.0 - 10.0; // -10.0 <-> 10.0
            }
        
            rkParameters<double> example(tripTime, alpha, beta, zeta, testcoeff); 
        
            inputParameters[i].startParams = example;
        }
    }
    // If not a random start, read from file
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
            int row = mt_rand() % numStarts; // Choose a random row to get the parameters from

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


    Individual *survivors = new Individual[SURVIVOR_COUNT]; // stores the winners of the head-to-head competition
    int newInd = numThreads; // the whole population is new the first time through the loop

    // setup output of generation results over time onto a .csv file
    std::ofstream generationPerformanceBestExcel;
    generationPerformanceBestExcel.open("BestInGenerations.csv");
    // Set first row in the file be a header for the columns
    generationPerformanceBestExcel << "Gen #" << "," << "posDiff" << "," << "velDiff" << "," 
                               << "rFinal" << "," << "thetaFinal" << "," << "zFinal" << "," << "vrFinal" << "," << "vthetaFinal" << "," << "vzFinal" << ","
                               << "rInitial" << "," << "thetaInitial" << "," << "zInitial" << ","<< "vrInitial" << "," << "vthetaInitial" << "," << "vzInitial" << ","
                               << "alpha" << "," << "beta" << "," << "zeta" << "," << "anneal" << "," << "tripTime" << "\n";


   std::ofstream generationPerformanceWorstExcel;
   generationPerformanceWorstExcel.open("WorstInGenerations.csv");
   // Set first row in the file be a header for the columns
    generationPerformanceWorstExcel << "Gen #" << "," << "posDiff" << "," << "velDiff" << "," 
                               << "rFinal" << "," << "thetaFinal" << "," << "zFinal" << "," << "vrFinal" << "," << "vthetaFinal" << "," << "vzFinal" << ","
                               << "rInitial" << "," << "thetaInitial" << "," << "zInitial" << ","<< "vrInitial" << "," << "vthetaInitial" << "," << "vzInitial" << ","
                               << "alpha" << "," << "beta" << "," << "zeta" << "," << "anneal" << "," << "tripTime" << "\n";

    // setup output of generation results over time onto a .bin file
    std::ofstream generationBestPerformanceBin("BestInGenerations.bin", std::ios::binary);
    std::ofstream generationWorstPerformanceBin("WorstInGenerations.bin", std::ios::binary);

    std::ofstream generationThrustBestExcel, generationThrustWorstExcel, generationThrustBestBin, generationThrustWorstBin;
    if (thrust.type != thruster<double>::NO_THRUST) {
        generationThrustBestExcel.open("BestThrustGens.csv");
        generationThrustBestExcel << "gen,gamma0,gamma1,gamma2,gamma3,gamma4,gamma5,gamma6,tau0,tau1,tau2,coast0,coast1,coast2,coast3,coast4";
        generationThrustWorstExcel.open("WorstThrustGens.csv");
        generationThrustWorstExcel << "gen,gamma0,gamma1,gamma2,gamma3,gamma4,gamma5,gamma6,tau0,tau1,tau2,coast0,coast1,coast2,coast3,coast4";

        generationThrustBestBin.open("BestThrustGens.bin", std::ios::binary);
        generationThrustWorstBin.open("WorstThrustGens.bin", std::ios::binary);
    }

    double generation = 0;    // A counter for number of generations calculated
    
    // A do-while loop that continues until it is determined that the pool of inputParameters has reached desired tolerance level for enough individuals (best_count)
    
    double currentDistance; // Contains value for how far away the best individual is from the tolerance value
    double tolerance = cConstants->pos_threshold; // Tolerance for what is an acceptable solution (currently just the position threshold which is furthest distance from the target allowed)
                                                // This could eventually take into account velocity too and become a more complex calculation
    double dRate = 1.0e-7;

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
                double alpha = PI * 2*((static_cast<double>(mt_rand()) / mt_rand.max()) - 0.5); // -PI <-> PI
                double beta  = PI * ((static_cast<double>(mt_rand()) / mt_rand.max())); // 0 <-> PI
                double zeta  = PI * ((static_cast<double>(mt_rand()) / mt_rand.max()) - 0.5); // -PI/2 <-> PI/2

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

                // Set to be a bad individual
                inputParameters[k].posDiff = 1.0;
                inputParameters[k].velDiff = 0.0;
             }

            // calculate its new cost function
            inputParameters[k].getCost(cConstants);
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
        
        Individual currentBest;
        if (static_cast<int>(generation) % cConstants->change_check == 0) { // Compare current best individual to that from CHANGE_CHECK many generations ago. If they are the same, change size of mutations
            currentBest = inputParameters[0];
          
            if ( !(changeInBest(previousBestPos, previousBestVel, currentBest, dRate)) ) { // previousBest starts at 0 to ensure changeInBest = true on generation 0
                currentAnneal = currentAnneal * cConstants->anneal_factor;
                std::cout << "\n new anneal: " << currentAnneal << std::endl;
                if(trunc(currentBest.posDiff/dRate)==0 || trunc(currentBest.velDiff/dRate)==0) { 
                    dRate = dRate/10; 
                    std::cout << "\nnew dRate: " << dRate << std::endl;
                }
            }
            previousBestPos = currentBest.posDiff;
            previousBestVel = currentBest.velDiff;
        }

        // Display a '.' to the terminal to show that a generation has been performed
        // This also serves to visually seperate the generation display on the terminal screen
        std::cout << '.';

        // Write the best and worst Individuals in every write_freq generations into the files to view progress over generations
        if (static_cast<int>(generation) % cConstants->write_freq == 0) {
            writeIndividualToFiles(generationPerformanceBestExcel, generationBestPerformanceBin, generation, inputParameters[0], new_anneal);
            writeIndividualToFiles(generationPerformanceWorstExcel, generationWorstPerformanceBin, generation, inputParameters[numThreads-1], new_anneal);

            if (thrust.type != thruster<double>::NO_THRUST) {
                writeThrustToFiles(generationThrustBestExcel, generationThrustBestBin, generation, inputParameters[0], cConstants);
                writeThrustToFiles(generationThrustWorstExcel, generationThrustWorstBin, generation, inputParameters[numThreads-1], cConstants);
            }
        }

        // Only call terminalDisplay every DISP_FREQ, not every single generation
        if ( static_cast<int>(generation) % cConstants->disp_freq == 0) {
            terminalDisplay(inputParameters[0], generation);
        }

        // Create a new generation and increment the generation counter
        newInd = crossover(survivors, inputParameters, SURVIVOR_COUNT, numThreads, new_anneal, cConstants, thrust);
        ++generation;
        
        // If the current distance is still higher than the tolerance we find acceptable, perform the loop again
    } while ( !allWithinTolerance(tolerance, inputParameters, generation, cConstants) );

    
    // output the best Individuals of the final generation, using writeTrajectoryToFile()
    // Files outputted allows plotting of solutions in matlab
    double *start = new double[OPTIM_VARS];
    double cost = 0;

    // Output to excel
    double annealPlacement = 0; //setting anneal to be a placeholder value that has no real meaning as there will be no next generation for anneal to impact
    
    // Write the final best and worst performing individuals to their respective files
    writeIndividualToFiles(generationPerformanceBestExcel, generationBestPerformanceBin, generation, inputParameters[0], annealPlacement);
    writeIndividualToFiles(generationPerformanceWorstExcel, generationWorstPerformanceBin, generation, inputParameters[numThreads-1], annealPlacement);

    // Write the best individuals with best_count in total outputted in seperate binary files
    for (int i = 0; i < cConstants->best_count; i++) {
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
        writeTrajectoryToFile(start, cost, i + 1, thrust, cConstants);
    }

    // Close the performance files now that the algorithm is finished
    generationPerformanceBestExcel.close();
    generationBestPerformanceBin.close();
    generationPerformanceWorstExcel.close();
    generationWorstPerformanceBin.close();

    if (thrust.type != thruster<double>::NO_THRUST) {
        generationThrustBestExcel.close();
        generationThrustWorstExcel.close();
        generationThrustBestBin.close();
        generationThrustWorstBin.close();
    }

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
    std::cout << "- Device name: " << prop.name << std::endl << std::endl;
    cudaSetDevice(0);
    
    cudaConstants const * cConstants = new cudaConstants("../Config_Constants/genetic.config"); // Declare the genetic constants used, with file path being used

    // cudaConstants const * compareConstants = new cudaConstants("../Config_Constants/genetic.config"); // Declaring a comparison cudaConstants to verify no changes to values after optimize
    // Display contents of cConstants resulting from reading the file
    std::cout << *cConstants << std::endl;

    //if ( !(sameConstants(*cConstants, *compareConstants))) {
    //   std::cout << "\nERROR - cConstants and compareConstants not the same at begginning!\n";
    //}

    // pre-calculate a table of Earth's position within possible mission time range
    //----------------------------------------------------------------
    // Define variables to be passed into EarthInfo
    double startTime = 15778800; // 0.5 year (s)
    double endTime = 78894000; // 2.5 years (s)
    double timeRes = 3600; // (s) position of earth is calculated for every hour

    launchCon = new EarthInfo(startTime, endTime, timeRes, cConstants); // a global variable to hold Earth's position over time

    //----------------------------------------------------------------
    // Define the number of threads/individuals that will be used in optimize
    int blockThreads = 32;
    int numThreads = 2880; // the number of cores on a Tesla k40
    //int numThreads = 1920; // 384 cores on K620 * 5 = 1920

    //std::ofstream efficiencyGraph; // for viewing how many runge-kuttas ran per second for each combination of threads per block and total threads 
    //efficiencyGraph.open("efficiencyGraph.csv");
    std::cout << std::endl << "running optimize() with " << blockThreads << " threads per block and " << numThreads << " total threads" << std::endl;
    
    thruster<double> thrust(cConstants);

    optimize(numThreads, blockThreads, cConstants, thrust);

    //efficiencyGraph << blockThreads << "," << numThreads << "," << calcPerS  << "\n";
    //efficiencyGraph.close();
    
    //if ( !(sameConstants(*cConstants, *compareConstants))) {
    //    std::cout << "\nERROR - cConstants had changed after optimize()!\n";
    //}


    delete launchCon;
    delete cConstants;
    // delete compareConstants;
    
    return 0;
}