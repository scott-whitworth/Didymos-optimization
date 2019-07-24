#include "geneitAlgorithm.h"
#include "individuals.h"

void selectWinners(individuals* pool, int selectionSize, individuals* survivors, individuals* losers)
{
    for(int i = 0; i<(selectionSize+1)/2; i++)
    {
        if (greater(pool[i],pool[i+1]))
        {
            survivors[i] = pool[i];
            losers[i]=pool[i+1];
        }
        else
        {
            survivors[i] = pool[i+1];
            losers[i]=pool[i];  
        }
    }
    
}