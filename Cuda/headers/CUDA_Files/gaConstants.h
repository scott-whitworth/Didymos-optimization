#ifndef gaConstants_h
#define gaConstants_h

// genetic algorithm constraints

#define generationsNum 17500 // total number of generations

#define SURVIVOR_COUNT 250 //360 // number of individuals to use for crossover each generation--MUST BE DIVISIBLE BY 2 TO PAIR OFF FOR CROSSOVER

#define MUTATION_RATE 0.15 // fraction of new offspring to mutate

#define DOUBLE_MUTATION_RATE 0.2 // fraction of mutations to mutate two genes instead of just one

#define TRIPLE_MUTATION_RATE 0.05 // fraction of mutations to mutate three genes instead of just one or two

#define ANNEAL_MAX .01; // max amount to mutate by as a fraction of the initial random range of a parameter

#define ANNEAL_MIN .001; // max amount to mutate by as a fraction of the initial random range of a parameter

#endif