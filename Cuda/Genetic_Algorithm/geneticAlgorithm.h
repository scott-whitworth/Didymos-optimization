#ifndef GENETICALGORITHM_H
#define GENETICALGORITHM_H

#include "individuals.h"


// takes Individuals from the pool to compete in a head-to-head competition
// input: 
    // pool: the overall population (should be in randomized order)
    // selectionSize: the number of survivors to return
// output: 
    //survivors: the winners of the competition
void selectWinners(Individual* pool, int selectionSize, Individual* survivors);

#include "geneticAlgorithm.cpp"
#endif