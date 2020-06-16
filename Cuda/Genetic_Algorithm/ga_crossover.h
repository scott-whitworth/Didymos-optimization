#ifndef GA_CROSSOVER_h
#define GA_CROSSOVER_h

// Collection of functions to perform crossover operations on rkParameters
// A crossover mask is an array of elements equal to the number of elements in the rkParameters list
// The mask defines what elements come from partent n and parent m
#include <random>

#include "../Runge_Kutta/rkParameters.h"
#include "individuals.h"
#include "thruster.h"

// creates a new Individual from two parents
// mask determines which parent each gene is inherited from
rkParameters<double> generateNewIndividual(const rkParameters<double> & p1, const rkParameters<double> & p2, const int mask[], thruster<double> thrust);

// generates new Individuals from the survivors(winners of competition) and replaces the worst Individuals in the pool(population) with these new ones
// returns the number of new Individuals put into the pool
int crossover(Individual *survivors, Individual *pool, int survivorSize, int poolSize, double annealing, thruster<double> thrust);

#include "ga_crossover.cpp"
#endif

