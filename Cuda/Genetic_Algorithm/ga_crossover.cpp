#include "../Runge_Kutta/rkParameters.h"
#include "../Config_Constants/config.h"
#include "ga_crossover.h"
#include <iostream>
#include <chrono>

#define SECONDS_IN_YEAR 365.25*24*3600 // Used with getRand for triptime mutation scale

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
// Output: pool is unchanged, survivors contains an array of size selectionSize of individuals that were quasi-randomly chosen
void selectSurvivors(Individual* pool, int selectionSize, Individual* survivors) {
    for(int i = 0; i < selectionSize; i++) {
        // While the array is a shuffled, when selecting a survivor make a neighbor comparison to choose the one with a lower cost (at least somewhat better choice)
        if ( pool[2*i] < pool[(2*i)+1] ) {
            survivors[i] = pool[2*i];
        }
        else {
            survivors[i] = pool[(2*i)+1];
        }
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
// Output: Returns rkParameter object that is new individual
rkParameters<double> generateNewIndividual(const rkParameters<double> & p1, const rkParameters<double> & p2, const int * mask, thruster<double>& thrust) {
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
                newInd.coeff.coast[i - COAST_OFFSET] = p2.coeff.coast[i - COAST_OFFSET]/2.0 + p1.coeff.coast[i - COAST_ARRAY_SIZE]/2.0;
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

    return newInd;    
}


// In a given Individual's parameters, mutate one gene gauranteed. Randomly decide to mutate a second gene or third gene some times
// mutate a gene by adding or subtracting a small, random value on a parameter property
// Input: p1 - rkParameter that is taken to be the mutation base
//        rng - random number generator to use
//        annealing - a scalar value on the max random number when mutating
//        cConstants - holds properties to use such as mutation rates and mutation scales for specific parameter property types
//        thrust - Used to check if the thruster needs to be taken into account
// Output: Returns rkParameter object that is the mutated version of p1
rkParameters<double> mutate(const rkParameters<double> & p1, std::mt19937_64 & rng, double annealing, const cudaConstants* cConstants, thruster<double>& thrust) {
    rkParameters<double> newInd = p1;
    int genesToMutate = 1; // number of genes to mutate
    int mutateChance = rng() % 100;

    if (mutateChance < cConstants->triple_mutation_rate * 100) {
        genesToMutate = 3;
    }
    else if (mutateChance < cConstants->double_mutation_rate * 100) {
        genesToMutate = 2;
    }

    int mutatedGenes[3]; // index of genes to mutate
    mutatedGenes[0] = rng() % OPTIM_VARS;

    if (genesToMutate > 1) {
        do {
            mutatedGenes[1] = rng() % OPTIM_VARS;
        } while (mutatedGenes[1] == mutatedGenes[0]); // make sure that each mutated gene is unique
    }

    if (genesToMutate > 2) {
        do {
            mutatedGenes[2] = rng() % OPTIM_VARS;
        } while (mutatedGenes[2] == mutatedGenes[0] || mutatedGenes[2] == mutatedGenes[1]); // make sure that each mutated gene is unique
    }

    for (int i = 0; i < genesToMutate; i++) {
        int mutatedValue = mutatedGenes[i]; // the gene to mutate
        // alter thrust coefficients only when using a thruster
        if (thrust.type != thruster<double>::NO_THRUST) {
            //check coeff
            if ( (mutatedValue >= GAMMA_OFFSET) && (mutatedValue <= (GAMMA_OFFSET + GAMMA_ARRAY_SIZE-1)) ) { // Gamma value
                newInd.coeff.gamma[mutatedValue] += getRand(cConstants->gamma_mutate_scale * annealing, rng);
            }
            else if ( (mutatedValue >= TAU_OFFSET) && (mutatedValue <= (TAU_OFFSET + TAU_ARRAY_SIZE-1))) { // Tau value
                newInd.coeff.tau[mutatedValue-TAU_OFFSET] += getRand(cConstants->tau_mutate_scale * annealing, rng);
            }
            else if (mutatedValue >= COAST_OFFSET && mutatedValue <= (COAST_OFFSET + COAST_ARRAY_SIZE-1)) { // Coast value
                newInd.coeff.coast[mutatedValue-COAST_OFFSET] += getRand(cConstants->coast_mutate_scale * annealing, rng);
            }
        }
        else if (mutatedValue == TRIPTIME_OFFSET) { // Time final
            newInd.tripTime += SECONDS_IN_YEAR*getRand(cConstants->triptime_mutate_scale * annealing, rng);
        }
        else if (mutatedValue == ZETA_OFFSET) { // Zeta
            newInd.zeta += getRand(cConstants->zeta_mutate_scale * annealing, rng);
        }
        else if (mutatedValue == BETA_OFFSET) { // Beta
            newInd.beta += getRand(cConstants->beta_mutate_scale * annealing, rng);
            // A check to ensure beta remains in value range 0 to pi
            if (newInd.beta < 0) {
                newInd.beta = 0;
            }
            else if (newInd.beta > M_PI) {
                newInd.beta = M_PI;
            }
        }
        else if (mutatedValue == ALPHA_OFFSET) { // Alpha
            newInd.alpha += getRand(cConstants->alpha_mutate_scale * annealing, rng);
        }
    }
    return newInd;
}


// Create a new individual, using two parents with a mask and also possible mutation occurring
// Input: pool - pointer array to Individuals that is where the new individual is stored
//        survivors - pointer array to Individuals to access the two parents from
//        mask - pointer array of maskValues used to decide on which property from which parent is acquired (or average of the two)
//        parent1Index, parent2Index - integer values for where in the survivor array to get the parents from
//        newIndiviudalIndex - integer value for where in the pool to store the new individual
//        annealing - If mutation is occurring, passed into mutate function
//        rng - Random number generator that is used to determine if mutation occurs and is passed into mutate function if it is occurring
//        cConstants - holds properties that is passed into mutate, also contains mutation_rate value that is used to determine if mutation will occur
// Output: pool[newIndiviudalIndex] contains a newly generated individual that is combination of survivor[parent1Index] and survivor[parent2Index] with possibly slight value changes in 1,2,3 variables
void mutateNewIndividual(Individual *pool, Individual *survivors, int * mask, int parent1Index, int parent2Index, int newIndividualIndex, double annealing, std::mt19937_64 & rng, const cudaConstants* cConstants, thruster<double>& thrust) {
    // Generate a new individual, using the parents and mask
    pool[newIndividualIndex] = Individual();
    pool[newIndividualIndex].startParams = generateNewIndividual(survivors[parent1Index].startParams, survivors[parent2Index].startParams, mask, thrust);

    // With percent chance of occurring, perform a mutation on the new individual's starting params
    if (rng() % 100 < cConstants->mutation_rate * 100) {
        pool[newIndividualIndex].startParams = mutate(pool[newIndividualIndex].startParams, rng, annealing, cConstants, thrust);
    }
}

// Method that creates a pair of new Individuals from a pair of other individuals and a mask
// Input: pool - pointer array to Individuals that is where the new pair of individuals are stored
//        survivors - pointer array to Individuals to access the two parents from
//        mask - pointer array of maskValues used to decide on which property from which parent is acquired (or average of the two)
//        newIndCount - value that tracks number of newly created indiviudals in the pool so far in the newGeneration process, also impacts where to put the new individuals in the pool
//        int parentsIndex - value for determing where the pair of parent survivors are selected (parent 1 is at parentsIndex, parent 2 is offset by +1)
//        annealing - double variable passed onto mutateNewIndividual
//        poolSize - length of the pool array
//        rng - random number generator passed on to mutateNewIndividual
//        cConstants - passed on to mutateNewIndividual
// Output: pool contains two newly created individuals at (poolSize - 1 - newIndCount) and (poolSize - 2 - newIndCount)
//         mask is flipped in polarity (refer to flipMask method)
//         newIndCount is incremented by +2
void generateChildrenPair(Individual *pool, Individual *survivors, int * mask, int& newIndCount, int parentsIndex, double annealing, int poolSize, std::mt19937_64 & rng, const cudaConstants* cConstants, thruster<double>& thrust) { 
    // Determine where the parents and the new individual being created are located in the pool
    int parent1Index = parentsIndex;
    int parent2Index = parentsIndex + 1;

    int newIndividualIndex = poolSize - 1 - newIndCount;    // The new indiviudal is located at the end of the pool, up the number of new individuals already created
    // Generate new offspring with mask
    mutateNewIndividual(pool, survivors, mask, parent1Index, parent2Index, newIndividualIndex, annealing, rng, cConstants, thrust);
    newIndCount++;

    // Get the opposite offspring from the mask by flipping the mask
    newIndividualIndex--; // Decrement newIndividualIndex value to access where the next individual must be as newIndCount has increased
    flipMask(mask);
    mutateNewIndividual(pool, survivors, mask, parent1Index, parent2Index, newIndividualIndex, annealing, rng, cConstants, thrust);
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
int newGeneration(Individual *survivors, Individual *pool, int survivorSize, int poolSize, double annealing, const cudaConstants* cConstants, thruster<double>& thrust) {
    std::mt19937_64 rng(cConstants->time_seed);

    int * mask = new int[OPTIM_VARS];

    int newIndCount = 0; // Number of new individuals created so far (initially none), used in navigating through the pool when creating new individuals and returned at end of function
    int numPairs = survivorSize / 2; // Value for how many pairs to use and produce in each loop (as one iteration through a loop produces a new pair)

    // Generate two offspring through each crossover method, total is 4 * survivorSize offspring in pool
    // Every loop needs to reset the mask as it is flipped from generateChildrenPair, for ones that use randomization it also keeps from having same mask for all new pairs

    // Loop for wholeRandom mask
    for (int i = 0; i < numPairs; i++) {
        crossOver_wholeRandom(mask, rng);
        generateChildrenPair(pool, survivors, mask, newIndCount, 2*i, annealing, poolSize, rng, cConstants, thrust);
    }
    // Loop for averaging mask
    for (int i = 0; i < numPairs; i++) {
        crossOver_average(mask);
        generateChildrenPair(pool, survivors, mask, newIndCount, 2*i, annealing, poolSize, rng, cConstants, thrust);
    }
    // 2 loops for bundleVars mask, two seperate loops resulting from carry over of past code, also will allow easier changes in masks used (right now just using two bundleVars instead of two different ones)
    for (int i = 0; i < numPairs; i++) {
        crossOver_bundleVars(mask, rng);
        generateChildrenPair(pool, survivors, mask, newIndCount, 2*i, annealing, poolSize, rng, cConstants, thrust);
    }
    for (int i = 0; i < numPairs; i++) {
        crossOver_bundleVars(mask, rng);
        generateChildrenPair(pool, survivors, mask, newIndCount, 2*i, annealing, poolSize, rng, cConstants, thrust);
    }

    delete [] mask;
    return newIndCount;
}
