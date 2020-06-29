#include "geneticAlgorithm.h"

void selectWinners(Individual* pool, int selectionSize, Individual* survivors) {
    for(int i = 0; i < selectionSize; i++) {
        if ( pool[2*i] < pool[(2*i)+1] ) {
            survivors[i] = pool[2*i];
        }
        else {
            survivors[i] = pool[(2*i)+1];
        }
    }   
}

double posCost(Individual* pool, int size) {
    return abs(pool[0].posDiff - pool[size-1].posDiff);
}

double velCost(Individual* pool, int size) {
    return abs(pool[0].velDiff - pool[size-1].velDiff);
}

bool converge(Individual* pool, cudaConstants* gConstants) {
    return (posConverge(pool, gConstants) && velConverge(pool, gConstants));
}

bool posConverge(Individual* pool, cudaConstants* gConstants) {
    return pool[0].posDiff < gConstants->v_impact;
}

bool velConverge(Individual* pool, cudaConstants* gConstants) {
    return pool[0].velDiff > gConstants->speed_threshold;
}