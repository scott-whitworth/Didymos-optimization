#ifndef GA_CROSSOVER_h
#define GA_CROSSOVER_h

// Collection of functions to perform crossover operations on rkParameters
// A crossover mask is an array of elements equal to the number of elements in the rkParameters list
// The mask defines what elements come from partent n and parent m
// [ {elements}     {coefficients}          {other parameters} ]
// [   0-5          6-14,15-19,20-24, 25,    26, 27                                             ]

//TODO: There must be a better way than this:
// Elements 6 + Coefficients: GAMMA_ARRAY_SIZE 9 + TAU_ARRAY_SIZE 5 + COAST_ARRAY_SIZE 5 + 1 + 2 (Wetmass + timeFinal)
#define RK_SIZE 28

#include "rkParameters.h"

rkParameters<double> generateNewIndividual(const rkParameters<double> & p1, const rkParameters<double> & p2, const int mask[]);

#endif

