#ifndef gaConstants_h
#define gaConstants_h

// genetic algorithm constraints

#define generationsNum 20000 // total number of generations

//#define SURVIVOR_COUNT 240 // number of individuals to use for crossover each generation--MUST BE DIVISIBLE BY 2 TO PAIR OFF FOR CROSSOVER
// 240 (survivors) / 2 (parents per pair) * 8 (offspring per pair) = 960 = half of 1920 --for k620 GPU

#define SURVIVOR_COUNT 360 // number of individuals to use for crossover each generation--MUST BE DIVISIBLE BY 2 TO PAIR OFF FOR CROSSOVER
// 360 (survivors) / 2 (parents per pair) * 8 (offspring per pair) = 1440 = half of 2880 --for k40 GPU

#define MUTATION_RATE 0.15 // fraction of new offspring to mutate

#define DOUBLE_MUTATION_RATE 0.2 // fraction of mutations to mutate two genes instead of just one

#define TRIPLE_MUTATION_RATE 0.05 // fraction of mutations to mutate three genes instead of just one or two

#define ANNEAL_MAX .01 // max amount to mutate by as a fraction of the initial random range of a parameter

#define ANNEAL_MIN .001 // max amount to mutate by as a fraction of the initial random range of a parameter

//#define POSITION_THRESH 1.0E-8 // threshold for determining weights of position and velocity in comparisons of runs in Individual::greaterInd() NOTE THIS VALUE HAS NOT BEEN DECIDED
#define POSITION_THRESH 1.0e-25 // to ignore velDiff entirely

#endif