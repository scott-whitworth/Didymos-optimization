#ifndef gaConstants_h
#define gaConstants_h

// genetic algorithm constraints

//#define generationsNum 150000 // total number of generationsIR OFF FOR CROSSOVER
// 360 (survivors) / 2 (parents per pair) * 8 (offspring per pair) = 1440 = half of 2880 --for k40 GPU
/*
#define MUTATION_RATE 0.15 // fraction of new offspring to mutate

#define DOUBLE_MUTATION_RATE 0.2 // fraction of mutations to mutate two genes instead of just one

#define TRIPLE_MUTATION_RATE 0.05 // fraction of mutations to mutate three genes instead of just one or two

//#define ANNEAL_MAX .001 // max amount to mutate by as a fraction of the initial random range of a parameter

//#define ANNEAL_MIN .0001 // max amount to mutate by as a fraction of the initial random range of a parameter

#define ANNEAL_INITIAL 0.001 // initial value for annealing, meant to replace the previously used calculation involving ANNEAL_MIN and ANNEAL_MAX with something more simple

// threshold for how close the spacecraft must be to the asteriod at end of its trajectory
#define POSITION_THRESH 1.0e-6 //

#define ANNEAL_FACTOR 0.25 // factor by which annealing is changed when there is no change in the best individual over 100 generations

#define WRITE_FREQ 100 // Determine how many generations between writing the progress of the best individual onto a .csv and .bin file

#define DISP_FREQ 200 // Determine how many generations between outputing contents onto the terminal screen

#define BEST_COUNT 10 // Number of individuals that will be recorded when the algorithm is finished

#define CHANGE_CHECK 100 // How often it checks for if the best individual has changed, used in the basis of Jeremy's method of anneal value dependent on if there was no change
*/
#endif