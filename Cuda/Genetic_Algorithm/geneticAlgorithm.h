#ifndef GENETICALGORITHM_H
#define GENETICALGORITHM_H

#include "individuals.h"
#include "../Config_Constants/config.h"

// takes Individuals from the pool to compete in a head-to-head competition
// input: 
    // pool: the overall population (should be in randomized order)
    // selectionSize: the number of survivors to return
// output: 
    //survivors: the winners of the competition
void selectWinners(Individual* pool, int selectionSize, Individual* survivors);

// Calculates the range of positions from the asteroid within the pool
// Input:
    // pool: the entire generation (must already be sorted)
    // size: the number of Individuals in the pool
// Output:
    // costRange: the difference in cost function between the least and greatest Individuals
double posCost(Individual* pool, int size);

// Calculates the range of velocities from the asteroid within the pool
// Input:
    // pool: the entire generation (must already be sorted)
    // size: the number of Individuals in the pool
// Output:
    // costRange: the difference in cost function between the least and greatest Individuals
double velCost(Individual* pool, int size);

// Determines whether the pool meets requirements for convergence
// Input:
    // pool: the entire generation (must already be sorted)
    // size: the number of Individuals in the pool
// Output:
    // convgFlag: Boolean indicating whether the pool has converged
bool converge(Individual* pool, const cudaConstants* gConstants);

// Determines whether the pools' posDiffs have converged
// Input:
    // pool: the entire generation (must already be sorted)
// Output:
    // convgFlag: Boolean indicating whether posDiffs have converged
bool posConverge(Individual* pool, const cudaConstants* gConstants);

// Determines whether the pool's velDiffs have converged
// Input:
    // pool: the entire generation (must already be sorted)
// Output:
    // convgFlag: Boolean indicating whether velDiffs have converged
bool velConverge(Individual* pool, const cudaConstants* gConstants);

#include "geneticAlgorithm.cpp"
#endif