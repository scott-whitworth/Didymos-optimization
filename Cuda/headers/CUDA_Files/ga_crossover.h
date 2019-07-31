#ifndef GA_CROSSOVER_h
#define GA_CROSSOVER_h

// Collection of functions to perform crossover operations on rkParameters
// A crossover mask is an array of elements equal to the number of elements in the rkParameters list
// The mask defines what elements come from partent n and parent m
// [ {elements}     {coefficients}          {other parameters} ]
// [   0-5          6-14,15-19,20-24, 25,    26, 27                                             ]

#include "rkParameters.h"
#include "individuals.h"
#include <random>

rkParameters<double> generateNewIndividual(const rkParameters<double> & p1, const rkParameters<double> & p2, const int mask[]);

// uses the different crossover methods to get offspring from a set of survivors AND replaces the worst Individuals in the pool with the new offspring
// the number of Indiviudals to replace may change, but that is all handled inside this function
// INPUT: survivors: the Individuals to crossover
//        pool: the total collection of current Individuals
//        survivorSize: the number of Individuals in survivors 
//        poolSize: the number of Individuals in pool
// RETURNS: the number of new individuals put in the pool
int crossover(Individual *survivors, Individual *pool, int survivorSize, int poolSize);

rkParameters<double> mutate(const rkParameters<double> & p1, mt19937_64 & rng);

#include "ga_crossover.cu"
#endif

