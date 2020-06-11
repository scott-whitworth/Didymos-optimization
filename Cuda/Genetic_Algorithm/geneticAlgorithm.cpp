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

bool converge(Individual* pool, size) {
    return posConverge(pool, size) && velConverge(pool, size);
}

bool posConverge(Individual* pool, int size) {
    if (pool[size-1].posDiff > IMPACT_THRESH) {
        if (posCost(pool, size)/pool[0].posDiff < CONVG_TOL) {
            return true;
        }
    }
    else return false;
}

bool velConverge(Individual* pool, int size) {
    if (pool[size-1].velDiff > SPEED_THRESH) {
        if (velCost(pool, size)/pool[0].velDiff < CONVG_TOL) {
            return true;
        }
    }
    else return false;
}