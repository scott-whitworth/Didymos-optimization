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

double posCost(Individual* pool) {
    return abs(pool[0].posDiff - pool[pool.size()].posDiff);
}

double velCost(Individual* pool) {
    return abs(pool[0].velDiff - pool[pool.size()].velDiff);
}

bool converge(Individual* pool) {
    return posConverge(pool) && velConverge(pool);
}

bool posConverge(Individual* pool) {
    if (pool[pool.size()-1].posDiff > IMPACT_THRESH) {
        if (abs(pool[0].posDiff - pool[pool.size()-1].posDiff)/pool[0].posDiff < CONVG_TOL) {
            return true;
        }
    }
    else return false;
}

bool velConverge(Individual* pool) {
    if (pool[pool.size()-1].velDiff > SPEED_THRESH) {
        if (abs(pool[0].velDiff - pool[pool.size()-1].velDiff)/pool[0].velDiff < CONVG_TOL) {
            return true;
        }
    }
    else return false;
}