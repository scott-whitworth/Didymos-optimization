// Collection of functions to perform crossover operations on rkParameters
// A crossover mask is an array of elements equal to the number of elements in the rkParameters list
// The mask defines what elements come from parent 1 and parent 2 or if it is an average of the two values

#include "../Runge_Kutta/rkParameters.h"
#include "../Config_Constants/config.h"
#include "ga_crossover.h"
#include <iostream>
#include <chrono>

#define SECONDS_IN_YEAR 365*24*3600 // Used with getRand for triptime mutation scale

// Global enumeration for the different mask values instead of 1,2,3 for better readibility
enum maskValue {
    PARTNER1 = 1,
    PARTNER2,
    AVG,
};

// Creates a random bifurcation mask
// Randomly picks one index to be the start of the '2's from mask
// input: rng - a constructed mt19937_64 random number generator
// in/out: mask - all data will be overwritten
//              - Based on random index, first selection will be 1's, last selection will be 2's
//              - ex: [1, 1, 1, 1, 2, 2]
void crossOver_randHalf(int * mask, int size, std::mt19937_64 & rng) {
    int crossIndex = rng() % (size-1);
    //cout << "Random Index: " << crossIndex << endl;
    for (int i = 0; i < size; i++) {
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
// Each element in a mask is randomly set to either 1 or 2
// in/out : All data overwritten, set randomly
// input: rng - a constructed mt19937_64 random number generator
void crossOver_wholeRandom(int * mask, int size, std::mt19937_64 & rng) {
    for (int i = 0; i < size; i++ ) {
        if (rng() % 2) { //Coin flip, either 1/0
            mask[i] = PARTNER2;
        }
        else {
            mask[i] = PARTNER1;
        }
    }
    return;
}

// This crossover method randomly chooses between partners for each variable. Similar to crossOver_wholeRandom, but this keeps all parameters for each variable together
void crossOver_bundleVars(int * mask, std::mt19937_64 & rng, cudaConstants * cConstants) {
    int p_gamma = 1 + rng() % 2; // partners for each variable are randomly chosen between 1 and 2
    int p_tau = 1 + rng() % 2;
    int p_coast = 1 + rng() % 2;
    int p_triptime = 1 + rng() % 2;
    int p_alpha = 1 + rng() % 2;
    int p_beta = 1 + rng() % 2;
    int p_zeta = 1 + rng() % 2;

    // Gamma
    for (int i = cConstants->gamma_offset; i < (cConstants->gamma_offset + cConstants->gamma_size); i++) {
        mask[i] = p_gamma;
    }
    // Tau
    for (int i = cConstants->tau_offset; i < (cConstants->tau_offset + cConstants->tau_size); i++) {
        mask[i] = p_tau;
    }
    // Coast
    for (int i = cConstants->coast_offset; i < (cConstants->coast_offset + cConstants->coast_size); i++) {
        mask[i] = p_coast;
    }

    mask[cConstants->triptime_offset] = p_triptime;

    mask[cConstants->alpha_offset] = p_alpha;
    
    mask[cConstants->beta_offset] = p_beta;
    
    mask[cConstants->zeta_offset] = p_zeta;

    return;
}

// Sets the entire mask to be AVG
void crossOver_average(int * mask, int size) {
    for (int i = 0; i < size; i++) {
        mask[i] = AVG;
    }
    return;
}

// CURRENTLY NOT IN USE. THIS FUNCTION AND crossOver_tauPos WERE REPLACED BY crossOver_bundleVars
//Create a mask the flips just the gamma coefficients
// input mask: over writes all data with [1 ... 1, 2, 2, ... 2, 2, 1 ... 1]
// 2's correspond to Gamma coefficients
/*void crossOver_gammaPos(int mask[], cudaConstants * cConstants) {
    for (int i = 0; i < cConstants->optim_vars; i++) {
        if ( (i >= cConstants->gamma_offset) && (i < (cConstants->gamma_offset + cConstants->gamma_size) ) ) {
            mask[i] = 2;
        } 
        else {
            mask[i] = 1;
        }
    }
    return;
}*/

// CURRENTLY NOT IN USE. THIS FUNCTION AND crossOver_tauPos WERE REPLACED BY crossOver_bundleVars
//Create a mask the flips just the tau coefficients
// input mask: over writes all data with [1 ... 1, 2, 2, ... 2, 2, 1 ... 1]
// 2's correspond to tau coefficients
/*void crossOver_tauPos(int mask[], cudaConstants * cConstants) {
    for (int i = 0; i < OPTIM_VARS; i++) {
        if ( (i >= cConstants->tau_offset) && (i < (cConstants->tau_offset + cConstants->tau_size) ) ) {
            mask[i] = 2;
        }
        else {
            mask[i] = 1;
        }
    }
    return;
}*/

//Utility to flip the polarity of a mask
// input:  mask is an array of size OPTIM_VARS, input with either 1 or 2 as a mask
// output: each 1 in mask will be reassigned to be a 2, each 2 will be reassigned 1
void flipMask(int * mask, int size) {
    for (int i = 0; i < size; i++) {
        if (mask[i] == 1) {
            mask[i] = 2;
        }
        else {
            mask[i] = 1;
        }
    }
    return;
}

//Copy contents of maskIn into maskOut
void copyMask(int * maskIn, int * maskOut, int size) {
    for (int i = 0; i < size; i++) {
        maskOut[i] = maskIn[i];
    }
}

// Display mask contents onto the terminal
void printMask(int * mask, int size) {
    std::cout << "[";

    for(int i = 0; i < size; i++) {
        std::cout << mask[i];
        // If not the last item, need a comma to seperate between the items in the display
        if ( i < size-1 ) {
            std::cout <<", ";
        }
    }
    std::cout << "]";
}

rkParameters<double> generateNewIndividual(const rkParameters<double> & p1, const rkParameters<double> & p2, const int * mask, thruster<double>& thrust, cudaConstants * cConstants) {
    // First set the new individual to hold traits from parent 1, then go through the mask to determine if a parameter value is to be set to parent 2 or be an average of parent 1 and 2
    rkParameters<double> newInd = p1;

    // Only calculate values pertaining to thrust if a thruster is being used
    if (thrust.type != thruster<double>::NO_THRUST) {
        // Iterate through the gamma values
        for (int i = cConstants->coast_offset; i < (cConstants->gamma_offset + cConstants->gamma_size); i++) {
            if (mask[i] == PARTNER2) {
                newInd.coeff.gamma[i - cConstants->gamma_offset] = p2.coeff.gamma[i - cConstants->gamma_offset]; // If the mask is 2, use the value from parent 2
            }
            else if (mask[i] == AVG) {
                newInd.coeff.gamma[i - cConstants->gamma_offset] = p2.coeff.gamma[i - cConstants->gamma_offset]/2.0 + p1.coeff.gamma[i - cConstants->gamma_offset]/2.0;// If the mask is 3, average the values from both parents
            }
        }
        // Iterate through tau values
        for (int i = cConstants->tau_offset; i < (cConstants->tau_offset + cConstants->tau_size); i++) {
            if (mask[i] == PARTNER2) {
                newInd.coeff.tau[i - cConstants->tau_offset] = p2.coeff.tau[i - cConstants->tau_offset];
            }
            else if (mask[i] == AVG) {
                newInd.coeff.tau[i - cConstants->tau_offset] = p2.coeff.tau[i - cConstants->tau_offset]/2.0 + p1.coeff.tau[i - cConstants->tau_offset]/2.0;
            }
        }
        // Iterate through coasting values
        for (int i = cConstants->coast_offset; i < (cConstants->coast_offset + cConstants->coast_size); i++) {
            if (mask[i] == PARTNER2) {
                newInd.coeff.coast[i - cConstants->coast_offset] = p2.coeff.coast[i - cConstants->coast_offset];
            }
            else if (mask[i] == AVG) {
                newInd.coeff.coast[i - cConstants->coast_offset] = p2.coeff.coast[i - cConstants->coast_offset]/2.0 + p1.coeff.coast[i - cConstants->coast_size]/2.0;
            }
        }
    }
    // Go through other variables
    if (mask[cConstants->triptime_offset] == PARTNER2) { //tripTime
        newInd.tripTime = p2.tripTime;
    }
    else if (mask[cConstants->triptime_offset] == AVG) {
        newInd.tripTime = p2.tripTime/2.0 + p1.tripTime/2.0;
    }
    if (mask[cConstants->zeta_offset] == PARTNER2) { //zeta
        newInd.zeta = p2.zeta;
    }
    else if (mask[cConstants->zeta_offset] == AVG) {
        newInd.zeta = p1.zeta/2.0 + p2.zeta/2.0;
    }
    if (mask[cConstants->beta_offset] == PARTNER2) { //beta
        newInd.beta = p2.beta;
    }
    else if (mask[cConstants->beta_offset] == AVG) {
        newInd.beta = p2.beta/2.0 + p1.beta/2.0;
    }
    if (mask[cConstants->alpha_offset] == PARTNER2) { //alpha
        newInd.alpha = p2.alpha;
    }
    else if (mask[cConstants->alpha_offset] == AVG) {
        newInd.alpha = p1.alpha/2.0 + p2.alpha/2.0;
    }

    return newInd;    
}

// utility function for mutate() to get a random double with as high a resolution as possible
// INPUT:
//      max: the absolute value of the min and max(min = -max) of the range
double getRand(double max, std::mt19937_64 & rng) {
    return static_cast<double>(rng()) / rng.max() * max * 2.0 - max;
}


// in a given Individual's parameters, mutate one gene gauranteed. Randomly decide to mutate a second gene some times.
// mutate a gene by adding or subtracting a small, random value from a parameter
rkParameters<double> mutate(const rkParameters<double> & p1, std::mt19937_64 & rng, double annealing, cudaConstants* cConstants, thruster<double>& thrust) {
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

    mutatedGenes[0] = rng() % cConstants->optim_vars;

    if (genesToMutate > 1) {
        do {
            mutatedGenes[1] = rng() % cConstants->optim_vars;
        } while (mutatedGenes[1] == mutatedGenes[0]); // make sure that each mutated gene is unique
    }

    if (genesToMutate > 2) {
        do {
            mutatedGenes[2] = rng() % cConstants->optim_vars;
        } while (mutatedGenes[2] == mutatedGenes[0] || mutatedGenes[2] == mutatedGenes[1]); // make sure that each mutated gene is unique
    }

    for (int i = 0; i < genesToMutate; i++) {
        int mutatedValue = mutatedGenes[i]; // the gene to mutate
        // alter thrust coefficients only when using a thruster
        if (thrust.type != thruster<double>::NO_THRUST) {
            //check coeff
            if ( (mutatedValue >= cConstants->gamma_offset) && (mutatedValue <= (cConstants->gamma_offset + cConstants->gamma_size-1)) ) {
                newInd.coeff.gamma[mutatedValue] += getRand(cConstants->gamma_mutate_scale * annealing, rng);
            }
            else if ( (mutatedValue >= cConstants->tau_offset) && (mutatedValue <= (cConstants->tau_offset + cConstants->tau_size-1))) {//Tau (7-9)
                newInd.coeff.tau[mutatedValue-7] += getRand(cConstants->tau_mutate_scale * annealing, rng);
            }
            else if (mutatedValue >= cConstants->coast_offset && mutatedValue <= (cConstants->coast_offset + cConstants->coast_size -1)) {//coast (14-18)
                newInd.coeff.coast[mutatedValue-14] += getRand(cConstants->coast_mutate_scale * annealing, rng);
            }
        }
        if (mutatedValue == cConstants->triptime_offset) { //Time final
            newInd.tripTime += SECONDS_IN_YEAR*getRand(cConstants->triptime_mutate_scale * annealing, rng);
        }
        if (mutatedValue == cConstants->zeta_offset) { //zeta
            newInd.zeta += getRand(cConstants->zeta_mutate_scale * annealing, rng);
        }
        if (mutatedValue == cConstants->beta_offset) { //beta
            newInd.beta += getRand(cConstants->beta_mutate_scale * annealing, rng);
        }
        if (mutatedValue == cConstants->alpha_offset) { //alpha
            newInd.alpha += getRand(cConstants->alpha_mutate_scale * annealing, rng);
        }
    }
    

    return newInd;    
}


// if using this, consider that theta is the same as theta + 2 pi
// CURRENTLY NOT USED
//Creates a new rkParameter based on the average between p1 and p2
// input: p1 and p2 are valid rkParameters
// output: average of the two
/*rkParameters<double> generateNewIndividual_avg(const rkParameters<double> & p1, const rkParameters<double> & p2, thruster<double>& thrust) {
    rkParameters<double> newInd;
    
    // alter thrust coefficients only when using a thruster
    if (thrust.type != thruster<double>::NO_THRUST) {
        for (int i = 0; i < p1.coeff.gammaSize; i++) {
            newInd.coeff.gamma[i] = (p1.coeff.gamma[i]/2.0) + (p2.coeff.gamma[i]/2.0);
        }
        for (int i = 0; i < p1.coeff.tauSize; i++) {
            newInd.coeff.tau[i] = (p1.coeff.tau[i]/2.0) + (p2.coeff.tau[i]/2.0);
        }
        for (int i = 0; i < p1.coeff.coastSize; i++) {
            newInd.coeff.coast[i] = (p1.coeff.coast[i]/2.0) + (p2.coeff.coast[i]/2.0);
        }
    }

    newInd.alpha = (p1.alpha/2.0) + (p2.alpha/2.0);
    newInd.beta = (p1.beta/2.0) + (p2.beta/2.0);
    newInd.zeta = (p1.zeta/2.0) + (p2.zeta/2.0);
    newInd.tripTime = (p1.tripTime/2.0) + (p2.tripTime/2.0);

    return newInd;    
}*/

// WIP
// Uses generateNewIndividual to create new Individuals by crossing over properties with mask, followed by random chance for mutations
// Output is 2 new individual in pool
// Can only be used within crossover function in its current state because of int i input
void mutateNewIndividual(Individual *pool, Individual *survivors, int mask[], int index, int i, double annealing, int poolSize, std::mt19937_64 & rng, cudaConstants* cConstants, thruster<double>& thrust) {
    pool[poolSize - 1 - (2 * index)] = Individual(); // create a new Individual instead of overwriting values
    pool[poolSize - 1 - (2 * index)].startParams = generateNewIndividual(survivors[2*i].startParams, survivors[(2*i)+1].startParams, mask, thrust, cConstants);
    
    if (rng() % 100 < cConstants->mutation_rate * 100) { // a certain chance of mutation
        pool[poolSize - 1 - (2 * index)].startParams = mutate(pool[poolSize - 1 - (4 * index)].startParams, rng, annealing, cConstants, thrust);
    }

    flipMask(mask, cConstants->optim_vars); // get the opposite offspring
    pool[poolSize - 1 - (2 * index) - 1] = Individual();
    pool[poolSize - 1 - (2 * index) - 1].startParams = generateNewIndividual(survivors[2*i].startParams, survivors[(2*i)+1].startParams, mask, thrust, cConstants);
    
    if (rng() % 100 < cConstants->mutation_rate * 100) {
        pool[poolSize - 1 - (2 * index) - 1].startParams = mutate(pool[poolSize - 1 - (4 * index) - 1].startParams, rng, annealing, cConstants, thrust);
    }
}


int crossover(Individual *survivors, Individual *pool, int survivorSize, int poolSize, double annealing, cudaConstants* cConstants, thruster<double>& thrust) {
    std::mt19937_64 rng(cConstants->time_seed);

    int * mask = new int[cConstants->optim_vars];

    int index = 0;

    // Generate two offspring through each crossover method, total of 8 offspring per parent pair
    for (int i = 0; i < survivorSize / 2; i++) {
        crossOver_wholeRandom(mask, cConstants->optim_vars, rng);
        mutateNewIndividual(pool, survivors, mask, index, i, annealing, poolSize, rng, cConstants, thrust);
        index++;
    }

    for (int i = 0; i < survivorSize / 2; i++) {
        crossOver_average(mask, cConstants->optim_vars);
        mutateNewIndividual(pool, survivors, mask, index, i, annealing, poolSize, rng, cConstants, thrust);
        index++;
    }

    for (int i = 0; i < survivorSize / 2; i++) {
        crossOver_bundleVars(mask, rng, cConstants);
        mutateNewIndividual(pool, survivors, mask, index, i, annealing, poolSize, rng, cConstants, thrust);
        index++;
    }

    for(int i = 0; i < survivorSize / 2; i++){
        crossOver_bundleVars(mask, rng, cConstants);
        mutateNewIndividual(pool, survivors, mask, index, i, annealing, poolSize, rng, cConstants, thrust);
        index++;
    }
    delete [] mask;
    return index * 2;
}
