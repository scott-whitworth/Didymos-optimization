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