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

double posCost(Individual* pool, int size) {
    return abs(pool[0].posDiff - pool[size-1].posDiff);
}

double velCost(Individual* pool, int size) {
    return abs(pool[0].velDiff - pool[size-1].velDiff);
}

bool converge(Individual* pool) {
    return posConverge(pool) && velConverge(pool);
}

bool posConverge(Individual* pool) {
    return pool[0].posDiff < POSITION_THRESH;
}

bool velConverge(Individual* pool) {
    return pool[0].velDiff > SPEED_THRESH;
}