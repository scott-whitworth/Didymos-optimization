#include "../Runge_Kutta/rkParameters.h"
#include "../Config_Constants/config.h"
#include "../output.h"
#include "ga_crossover.h"
#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>

// Global enumeration for the different mask values instead of 1,2,3 for better readibility and clarity of value meaning
enum maskValue {
    PARTNER1 = 1,
    PARTNER2,
    AVG,
};

// Method of determing selection of survivors that will carry properties into the new individuals of the newGeneration
// Input: pool - a shuffled pointer array of individuals to choose from
//        selectionSize - integer number of how many survivors to choose out of the pool
//        survivors - pointer array of individuals to copy the selected individuals and store
// Output: pool is unchanged, survivors contains an array of size selectionSize of individuals that contains survivors that are 
void selectSurvivors(Individual * pool, int poolSize, int selectionSize, Individual* survivors) {
    // Sort the pool by positional difference and make half the selctions the best posDiff
    std::sort(pool, pool+poolSize, BetterPosDiff);
    for (int i = 0; i < selectionSize / 2; i++) {
        survivors[i] = pool[i];
    }

    // Sort the pool by positional difference
    std::sort(pool, pool+poolSize, BetterVelDiff);

    for (int i = 0; i < selectionSize / 2; i++) {
        survivors[i] = pool[i];
    }


    return;
}

// Creates a random bifurcation mask, currently not in use
// Randomly picks one index to be the start of the '2's from mask
// input: rng - a constructed mt19937_64 random number generator
// in/out: mask - all data will be overwritten
//              - Based on random index, first selection will be PARTNER1's, last selection will be PARTNER2's
//              - ex: [1, 1, 1, 1, 2, 2]
void crossOver_randHalf(int * mask, std::mt19937_64 & rng) {
    int crossIndex = rng() % (OPTIM_VARS-1);
    //cout << "Random Index: " << crossIndex << endl;
    for (int i = 0; i < OPTIM_VARS; i++) {
        if (i > crossIndex) {
            mask[i] = PARTNER2;
        }
        else {
            mask[i] = PARTNER1;
        }        
    }
    return;
}

// Sets the entire mask to be PARTNER1 for length OPTIM_VARS, allows a crossover where no mixing occurs, currently not in use
// Input: mask - pointer integer array of length OPTIM_VARS
// Output: mask is set to contain all PARTNER1 values
void crossOver_oneParent(int * mask) {
    for (int i = 0; i < OPTIM_VARS; i++) {
        mask[i] = PARTNER1;
    }
}

// Creates a random mask
// Each element in a mask is randomly set to either PARTNER1 or PARTNER2
// in/out : All data overwritten, set randomly
// input: rng - a constructed mt19937_64 random number generator
void crossOver_wholeRandom(int * mask, std::mt19937_64 & rng) {
    for (int i = 0; i < OPTIM_VARS; i++ ) {
        if (rng() % 2) { //Coin flip, either 1/0
            mask[i] = PARTNER2;
        }
        else {
            mask[i] = PARTNER1;
        }
    }
    return;
}

// This crossover method randomly chooses between partners for each variable. Similar to crossOver_wholeRandom, but this keeps all parameters types (gamma, tau, coast) grouped
// Input: int pointer array mask of size OPTIM_VARS, rng that is a constructed mt19937_64 random number generator used to derive random values between 1 and 2
// Output: mask will contain an array of values that is either 1 or 2 (equivalent to PARTNER1 or PARTNER2)
void crossOver_bundleVars(int * mask, std::mt19937_64 & rng) {
    int p_gamma = 1 + rng() % 2; // partners for each variable are randomly chosen between 1 and 2
    int p_tau = 1 + rng() % 2;
    int p_coast = 1 + rng() % 2;
    int p_triptime = 1 + rng() % 2;
    int p_alpha = 1 + rng() % 2;
    int p_beta = 1 + rng() % 2;
    int p_zeta = 1 + rng() % 2;

    // Gamma values
    for (int i = GAMMA_OFFSET; i < (GAMMA_OFFSET + GAMMA_ARRAY_SIZE); i++) {
        mask[i] = p_gamma;
    }
    // Tau values
    for (int i = TAU_OFFSET; i < (TAU_OFFSET + TAU_ARRAY_SIZE); i++) {
        mask[i] = p_tau;
    }
    // Coast values
    for (int i = COAST_OFFSET; i < (COAST_OFFSET + COAST_ARRAY_SIZE); i++) {
        mask[i] = p_coast;
    }

    mask[TRIPTIME_OFFSET] = p_triptime;
    mask[ALPHA_OFFSET] = p_alpha;
    mask[BETA_OFFSET] = p_beta;
    mask[ZETA_OFFSET] = p_zeta;

    return;
}

// Sets the entire mask to be AVG for length OPTIM_VARS
// Input: mask - pointer integer array of length OPTIM_VARS
// Output: mask is set to contain all AVG values
void crossOver_average(int * mask) {
    for (int i = 0; i < OPTIM_VARS; i++) {
        mask[i] = AVG;
    }
    return;
}


// Utility to flip the polarity of a mask
// Input:  mask is an array of size OPTIM_VARS, input based on maskValue enumerations as a mask
// Output: each PARTNER1 in mask will be reassigned to be a PARTNER2, each PARTNER2 will be reassigned PARTNER1, AVG will be unchanged
void flipMask(int * mask) {
    for (int i = 0; i < OPTIM_VARS; i++) {
        if (mask[i] == PARTNER1) {
            mask[i] = PARTNER2;
        }
        else if (mask[i] == PARTNER2) {
            mask[i] = PARTNER1;
        }
        // If mask[i] is neither partner1 nor partner2 (must be avg then), leave it be
    }
    return;
}

// Copy contents of maskIn into maskOut of size OPTIM_VARS
// Input: maskIn an array that is the original mask, maskOut is an array that will have contents copied to
// Output: maskIn remains unchanged, maskOut will have same contents as maskIn
void copyMask(int * maskIn, int * maskOut) {
    for (int i = 0; i < OPTIM_VARS; i++) {
        maskOut[i] = maskIn[i];
    }
}

// Utility function for mutate() to get a random double with high resolution
// Input: max - the absolute value of the min and max(min = -max) of the range
// Output: A double value that is between -max and +max
double getRand(double max, std::mt19937_64 & rng) {
    return static_cast<double>(rng()) / rng.max() * max * 2.0 - max;
}

// Creates a new rkParameters individual by combining properties of two parent Individuals using a mask to determine which
// Input: two rkParameter individuals (p1 and p2) - used in creating the new individual
//        mask - Contains maskValue values and length of OPTIM_VARS, determines how the two parent properties are merged into creating the new individual
//        thrust - Determine if the thruster properties must be carried over
//        cConstants - passed to mutate()
//        annealing - passed to mutate()
//        rng - passed to mutate()
//        generation - passed to mutate()
// Output: Returns rkParameter object that is new individual
rkParameters<double> generateNewIndividual(const rkParameters<double> & p1, const rkParameters<double> & p2, const int * mask, thruster<double>& thrust, const cudaConstants * cConstants, double annealing, std::mt19937_64 & rng, double generation) {
    // First set the new individual to hold traits from parent 1, then go through the mask to determine if a parameter value is to be set to parent 2 or be an average of parent 1 and 2
    rkParameters<double> newInd = p1;

    // Only calculate values pertaining to thrust if a thruster is being used
    if (thrust.type != thruster<double>::NO_THRUST) {
        // Iterate through the gamma values
        for (int i = GAMMA_OFFSET; i < (GAMMA_OFFSET + GAMMA_ARRAY_SIZE); i++) {
            if (mask[i] == PARTNER2) {
                newInd.coeff.gamma[i - GAMMA_OFFSET] = p2.coeff.gamma[i - GAMMA_OFFSET]; // If the mask is 2, use the value from parent 2
            }
            else if (mask[i] == AVG) {
                newInd.coeff.gamma[i - GAMMA_OFFSET] = p2.coeff.gamma[i - GAMMA_OFFSET]/2.0 + p1.coeff.gamma[i - GAMMA_OFFSET]/2.0;// If the mask is 3, average the values from both parents
            }
        }
        // Iterate through tau values
        for (int i = TAU_OFFSET; i < (TAU_OFFSET + TAU_ARRAY_SIZE); i++) {
            if (mask[i] == PARTNER2) {
                newInd.coeff.tau[i - TAU_OFFSET] = p2.coeff.tau[i - TAU_OFFSET];
            }
            else if (mask[i] == AVG) {
                newInd.coeff.tau[i - TAU_OFFSET] = p2.coeff.tau[i - TAU_OFFSET]/2.0 + p1.coeff.tau[i - TAU_OFFSET]/2.0;
            }
        }
        // Iterate through coasting values
        for (int i = COAST_OFFSET; i < (COAST_OFFSET + COAST_ARRAY_SIZE); i++) {
            if (mask[i] == PARTNER2) {
                newInd.coeff.coast[i - COAST_OFFSET] = p2.coeff.coast[i - COAST_OFFSET];
            }
            else if (mask[i] == AVG) {
                newInd.coeff.coast[i - COAST_OFFSET] = p2.coeff.coast[i - COAST_OFFSET]/2.0 + p1.coeff.coast[i - COAST_OFFSET]/2.0;
            }
        }
    }
    // Go through other variables
    if (mask[TRIPTIME_OFFSET] == PARTNER2) { //tripTime
        newInd.tripTime = p2.tripTime;
    }
    else if (mask[TRIPTIME_OFFSET] == AVG) {
        newInd.tripTime = p2.tripTime/2.0 + p1.tripTime/2.0;
    }
    if (mask[ZETA_OFFSET] == PARTNER2) { //zeta
        newInd.zeta = p2.zeta;
    }
    else if (mask[ZETA_OFFSET] == AVG) {
        newInd.zeta = p1.zeta/2.0 + p2.zeta/2.0;
    }
    if (mask[BETA_OFFSET] == PARTNER2) { //beta
        newInd.beta = p2.beta;
    }
    else if (mask[BETA_OFFSET] == AVG) {
        newInd.beta = p2.beta/2.0 + p1.beta/2.0;
    }
    if (mask[ALPHA_OFFSET] == PARTNER2) { //alpha
        newInd.alpha = p2.alpha;
    }
    else if (mask[ALPHA_OFFSET] == AVG) {
        newInd.alpha = p1.alpha/2.0 + p2.alpha/2.0;
    }
    // Now that newInd contains crossovered parameter values, call mutate onto it
    newInd = mutate(newInd, rng, annealing, cConstants, thrust, generation);

    return newInd;    
}

// Utility function to generate a boolean mask that determines which parameter value is mutating and how many based on mutation_rate iteratively
// input: rng - random number generating object used to randomly generate index values
//        mutateMask - pointer to a boolean array that is assumed length of OPTIM_VARS, sets genes being mutated to true and others to false
//        mutation_rate - a double value  less than 1 that is the chance a gene will be mutated, called iteratively to mutate more genes
// output: mutateMask contains false for genes that are not mutating, true for genes that are to be mutated
void mutateMask(std::mt19937_64 & rng, bool * mutateMask, double mutation_rate) {
    for (int i = 0; i < OPTIM_VARS; i++) {
        mutateMask[i] = false;
    }
    int geneCount = 0; // counter to make sure in the unlikely event that all genes are being mutated the code doesn't get stuck in infinite loop looking for a false spot in the array
    // Set a gene to mutate if a randomized values is less than mutation_rate, repeating everytime this is true
    while ((static_cast<double>(rng()) / rng.max()) < mutation_rate && geneCount < OPTIM_VARS) {
        bool geneSet = false; // boolean to flag if a gene was successfully selected to be set as true
        int index; // index value that will be assigned randomly to select a gene
        while (geneSet != true) {
            index = rng() % OPTIM_VARS;
            // If the randomly selected gene hasn't been already set to mutate, set it and flag the geneSet to true to break out of the loop
            if (mutateMask[index] == false) {
                mutateMask[index] = true;
                geneSet = true;
                geneCount++;
            }
        }

    }
}

// In a given Individual's parameters, generate a mutate mask using mutateMask() and then adjust parameters based on the mask, mutation of at least one gene is not guranteed
// mutate a gene by adding or subtracting a small, random value on a parameter property
// Input: p1 - rkParameter that is taken to be the mutation base
//        rng - random number generator to use
//        annealing - a scalar value on the max random number when mutating
//        cConstants - holds properties to use such as mutation rates and mutation scales for specific parameter property types
//        thrust - Used to check if the thruster needs to be taken into account
// Output: Returns rkParameter object that is the mutated version of p1
rkParameters<double> mutate(const rkParameters<double> & p1, std::mt19937_64 & rng, double annealing, const cudaConstants* cConstants, thruster<double>& thrust, double generation) {    
    rkParameters<double> newInd = p1; // initially set new individual to have all parameter values from parent 1

    // Declare and set a mutation_mask for which gene is being mutated
    bool * mutation_mask = new bool[OPTIM_VARS];
    mutateMask(rng, mutation_mask, cConstants->mutation_rate);

    // Declare a record that is to describe what genes are being changed and by how much to record into mutateFile
    double recordLog[OPTIM_VARS];

    // Iterate through the mutation_mask, mutating the corresponding gene if set to true
    for (int index = 0; index < OPTIM_VARS; index++) {
        if (mutation_mask[index] == true) {
            
            if ( (index >= GAMMA_OFFSET) && (index <= (GAMMA_OFFSET + GAMMA_ARRAY_SIZE-1)) ) { // Gamma value
                double randVar = getRand(cConstants->gamma_mutate_scale * annealing, rng);
                newInd.coeff.gamma[index-GAMMA_OFFSET] += randVar;
                recordLog[index] = randVar;
            }
            else if ( (index >= TAU_OFFSET) && (index <= (TAU_OFFSET + TAU_ARRAY_SIZE-1))) { // Tau value 
                double randVar = getRand(cConstants->tau_mutate_scale * annealing, rng);
                newInd.coeff.tau[index-TAU_OFFSET] += randVar;
                recordLog[index] = randVar;
            }
            else if (index >= COAST_OFFSET && index <= (COAST_OFFSET + COAST_ARRAY_SIZE-1)) { // Coast value
                double randVar = getRand(cConstants->coast_mutate_scale * annealing, rng);
                newInd.coeff.coast[index-COAST_OFFSET] += randVar;
                recordLog[index] = randVar;
            }
            else if (index == TRIPTIME_OFFSET) { // Time final
                double randVar = SECONDS_IN_YEAR * getRand(cConstants->triptime_mutate_scale * annealing, rng);
                newInd.tripTime += randVar;
                // bound checking to make sure the tripTime is set within the valid range of trip times
                if (newInd.tripTime < cConstants->triptime_min * SECONDS_IN_YEAR) {
                    newInd.tripTime = cConstants->triptime_min * SECONDS_IN_YEAR;
                }
                else if (newInd.tripTime > (cConstants->triptime_max) * SECONDS_IN_YEAR) {
                    newInd.tripTime = cConstants->triptime_max * SECONDS_IN_YEAR;
                }

                recordLog[index] = randVar;
            }
            else if (index == ZETA_OFFSET) { // Zeta
                double randVar = getRand(cConstants->zeta_mutate_scale * annealing, rng);
                newInd.zeta += randVar;
                recordLog[index] = randVar;
            }
            else if (index == BETA_OFFSET) { // Beta
                double randVar = getRand(cConstants->beta_mutate_scale * annealing, rng);
                newInd.beta = randVar;
                recordLog[index] = randVar;
    
                // A check to ensure beta remains in value range 0 to pi, doesn't update recordLog
                if (newInd.beta < 0) {
                    newInd.beta = 0;
                }
                else if (newInd.beta > M_PI) {
                    newInd.beta = M_PI;              
                }
            }
            else if (index == ALPHA_OFFSET) { // Alpha
                double randVar = getRand(cConstants->alpha_mutate_scale * annealing, rng);
                newInd.alpha += randVar;                
                recordLog[index] = randVar;
            }
        }
        else { // Record if the gene is not being mutated
            recordLog[index] = 0;
        }
    }

    // If in record mode, append the recordLog into the .csv file
    if (cConstants->record_mode == true) {
        int genesMutated = 0;
        for (int i = 0; i < OPTIM_VARS; i++) {
            if (mutation_mask[i] == true) {
                genesMutated++;
            }
        }
        recordMutateFile(cConstants, generation, annealing, genesMutated, recordLog);
    }
    delete [] mutation_mask;
    return newInd;
}


// Method that creates a pair of new Individuals from a pair of other individuals and a mask
// Input: pool - pointer array to Individuals that is where the new pair of individuals are stored
//        survivors - pointer array to Individuals to access the two parents from
//        mask - pointer array of maskValues used to decide on which property from which parent is acquired (or average of the two)
//        newIndCount - value that tracks number of newly created indiviudals in the pool so far in the newGeneration process, also impacts where to put the new individuals in the pool
//        int parentsIndex - value for determing where the pair of parent survivors are selected (parent 1 is at parentsIndex, parent 2 is offset by +1)
//        annealing - double variable passed onto generateNewIndividual
//        poolSize - length of the pool array
//        rng - random number generator passed on to generateNewIndividual
//        cConstants - passed on to generateNewIndividual
// Output: pool contains two newly created individuals at (poolSize - 1 - newIndCount) and (poolSize - 2 - newIndCount)
//         mask is flipped in polarity (refer to flipMask method)
//         newIndCount is incremented by +2
void generateChildrenPair(Individual *pool, Individual *survivors, int * mask, int& newIndCount, int parentsIndex, double annealing, int poolSize, std::mt19937_64 & rng, const cudaConstants* cConstants, thruster<double>& thrust, double generation) { 
    // Determine where the parents and the new individual being created are located in the pool
    int parent1Index = parentsIndex;
    int parent2Index = parentsIndex + 1;

    int newIndividualIndex = poolSize - 1 - newIndCount;    // The new indiviudal is located at the end of the pool, up the number of new individuals already created
    // Generate new offspring with mask
    pool[newIndividualIndex] = Individual(generateNewIndividual(survivors[parent1Index].startParams, survivors[parent2Index].startParams, mask, thrust, cConstants, annealing, rng, generation), cConstants);
    newIndCount++;

    // Get the opposite offspring from the mask by flipping the mask
    newIndividualIndex--; // Decrement newIndividualIndex value to access where the next individual must be as newIndCount has increased
    flipMask(mask);
    pool[newIndividualIndex] = Individual(generateNewIndividual(survivors[parent1Index].startParams, survivors[parent2Index].startParams, mask, thrust, cConstants, annealing, rng, generation), cConstants);
    newIndCount++;

    return;
}

// Creates the next pool to be used in the optimize function in opimization.cu
// Input: survivors - Individual pointer array of Individuals to be used in creating new individuals
//        pool - Sorted Individual pointer array that contains current generation
//        survivorSize - length of survivors array
//        poolSize - length of pool array
//        annealing - passed onto generateChildrenPair
//        cConstants - passed onto generateChildrenPair
//        thrust - passed onto generateChildrenPair
// Output: lower (survivorSize * 4) portion of pool is replaced with new individuals
//         Returns number of new individuals created (newIndCount)
int newGeneration(Individual *survivors, Individual *pool, int survivorSize, int poolSize, double annealing, const cudaConstants* cConstants, thruster<double>& thrust, std::mt19937_64 & rng, double generation) {

    int * mask = new int[OPTIM_VARS];
    int newIndCount = 0; // Number of new individuals created so far (initially none), used in navigating through the pool when creating new individuals and returned at end of function
    int numPairs = survivorSize / 2; // Value for how many pairs to use and produce in each loop (as one iteration through a loop produces a new pair)

    // Generate two offspring through each crossover method, total is 4 * survivorSize offspring in pool
    // Every loop needs to reset the mask as it is flipped from generateChildrenPair, for ones that use randomization it also keeps from having same mask for all new pairs

    // Loop for wholeRandom mask
    for (int i = 0; i < numPairs; i++) {
        crossOver_wholeRandom(mask, rng);
        generateChildrenPair(pool, survivors, mask, newIndCount, 2*i, annealing, poolSize, rng, cConstants, thrust, generation);
    }
    // Loop for averaging mask
    for (int i = 0; i < numPairs; i++) {
        crossOver_average(mask);
        generateChildrenPair(pool, survivors, mask, newIndCount, 2*i, annealing, poolSize, rng, cConstants, thrust, generation);
    }
    // 2 loops for bundleVars mask, two seperate loops resulting from carry over of past code, also will allow easier changes in masks used (right now just using two bundleVars instead of two different ones)
    for (int i = 0; i < numPairs; i++) {
        crossOver_bundleVars(mask, rng);
        generateChildrenPair(pool, survivors, mask, newIndCount, 2*i, annealing, poolSize, rng, cConstants, thrust, generation);
    }
    for (int i = 0; i < numPairs; i++) {
        crossOver_bundleVars(mask, rng);
        generateChildrenPair(pool, survivors, mask, newIndCount, 2*i, annealing, poolSize, rng, cConstants, thrust, generation);
    }
    delete [] mask;
    return newIndCount;
}
