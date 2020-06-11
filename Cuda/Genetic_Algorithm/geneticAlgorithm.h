#ifndef GENETICALGORITHM_H
#define GENETICALGORITHM_H

#include "individuals.h"
#include "gaConstants.h"


// takes Individuals from the pool to compete in a head-to-head competition
// input: 
    // pool: the overall population (should be in randomized order)
    // selectionSize: the number of survivors to return
// output: 
    //survivors: the winners of the competition
void selectWinners(Individual* pool, int selectionSize, Individual* survivors);

// Calculates the cost function range across the pool
// Input:
    // pool: the entire generation (must already be sorted)
    // size: the number of Individuals in the pool
// Output:
    // costRange: the difference in cost function between the least and greatest Individuals
double calcCost(Individual* pool, int size);

// Determines whether the pool meets requirements for convergence
// Input:
    // pool: the entire generation (must already be sorted)
    // size: the number of Individuals in the pool
// Output:
    // convgFlag: Boolean indicating whether the pool has converged
bool converge(Individual* pool, int size);

#include "geneticAlgorithm.cpp"
#endif