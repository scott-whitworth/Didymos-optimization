#ifndef gaConstants_h
#define gaConstants_h

// genetic algorithm constraints

#define SURVIVOR_COUNT 250 // number of individuals to use for crossover each generation--MUST BE DIVISIBLE BY 2 TO PAIR OFF FOR CROSSOVER

#define MUTATION_RATE 0.15 // percent of new offspring to mutate

#define DOUBLE_MUTATION_RATE 0.2 // percent chance for two genes to be mutated in an individual instead of just one

#define TRIPLE_MUTATION_RATE 0.05 // percent chance for three genes to be mutated in an individual instead of just one or two

#endif