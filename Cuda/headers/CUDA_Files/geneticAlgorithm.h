#ifndef GENETICALGORITHM_H
#define GENETICALGORITHM_H

#include "individuals.h"

// Assings in an array the winners and in a sepparte4 array the losers using the function from individuals greater to comparte
// input
// output
void selectWinners(Individual* pool, int selectionSize, Individual* survivors);

#include "geneticAlgorithm.cpp"
#endif