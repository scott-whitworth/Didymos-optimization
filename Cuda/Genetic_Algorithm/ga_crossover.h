#ifndef GA_CROSSOVER_h
#define GA_CROSSOVER_h

// Collection of functions to perform crossover operations on rkParameters
// A crossover mask is an array of elements equal to the number of elements in the rkParameters list
// The mask defines what elements come from partent n and parent m
#include <random>

#include "../Runge_Kutta/rkParameters.h"
#include "individuals.h"
#include "../Thrust_Files/thruster.h"
#include "../Config_Constants/config.h"

// takes Individuals from the pool to compete in a head-to-head competition
// input: 
    // pool: the overall population (should be in randomized order)
    // selectionSize: the number of survivors to return
// output: 
    //survivors: the winners of the competition
void selectWinners(Individual* pool, int selectionSize, Individual* survivors);

// creates a new Individual from two parents
// mask determines which parent each gene is inherited from
rkParameters<double> generateNewIndividual(const rkParameters<double> & p1, const rkParameters<double> & p2, const int mask[], thruster<double>& thrust);

// Has random chance of changing the parameters by a random quantity
rkParameters<double> mutate(const rkParameters<double> & p1, std::mt19937_64 & rng, double annealing, cudaConstants& gConstant, thruster<double>& thrust);

// Currently not in use, creates a new individual by averaging the parameter values between two parents
rkParameters<double> generateNewIndividual_avg(const rkParameters<double> & p1, const rkParameters<double> & p2);

// Uses mutate and generateNewIndividual to produce a new individual with chance of mutation and assigns them within the pool
void mutateNewIndividual(Individual *pool, Individual *survivors, int mask[], int index, int i, double annealing, int poolSize, std::mt19937_64 & rng, const cudaConstants* cConstants, thruster<double>& thrust);

// generates new Individuals from the survivors(winners of competition) and replaces the worst Individuals in the pool(population) with these new ones
// Uses mutateNewIndividual and the crossOver mask methods
// returns the number of new Individuals put into the pool
int newGeneration(Individual *survivors, Individual *pool, int survivorSize, int poolSize, double annealing, const cudaConstants* gConstant, thruster<double>& thrust);

void crossOver_randHalf(int mask[], std::mt19937_64 & rng);

void crossOver_wholeRandom(int mask[], std::mt19937_64 & rng);

void crossOver_gammaPos(int mask[]);

void crossOver_tauPos(int mask[]);

void flipMask(int mask[]);

void copyMask(int maskIn[], int maskOut[]);

void printMask(int mask[]);

#include "ga_crossover.cpp"
#endif

