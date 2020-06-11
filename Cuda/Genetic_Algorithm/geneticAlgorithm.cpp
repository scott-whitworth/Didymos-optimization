#include "geneticAlgorithm.h"

void selectWinners(Individual* pool, int selectionSize, Individual* survivors) {
    for(int i = 0; i < selectionSize; i++) {
        if (greaterInd(pool[2*i],pool[(2*i)+1])) {
            survivors[i] = pool[2*i];
        }
        else {
            survivors[i] = pool[(2*i)+1];
        }
    }   
}

double calcCost(Individual* pool, int size) {
    if (pool[size-1].posDiff >= IMPACT_THRESH) {
        return 1/AU;
    }
    else return abs(pool[0].velDiff - pool[size - 1].velDiff);
}