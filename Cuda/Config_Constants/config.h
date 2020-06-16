#ifndef GENETICCONFIG_h
#define GENETICCONFIG_h


// Structure that holds constant values related/used for the genetic algorithm that can be configured within genetic.config file
struct geneticConstants {
    double time_seed; // Seed used for randomization within optimize function, if it's set to -1 the seed is set to time(0) for genuine randomness
    double pos_threshold; // threshold for how close the spacecraft must be to the asteriod at end of its trajectory
    double anneal_factor; // factor by which annealing is changed when there is no change in the best individual over 100 generations

    int write_freq; // Determine how many generations between writing the progress of the best individual onto a .csv and .bin file
    int disp_freq; // Determine how many generations between outputing contents onto the terminal screen

    int best_count; // Number of individuals that needs to be within the acceptable condition before ending the algorithm, also how many of the top individuals are recorded
    int change_check; // How often it checks for if the best individual has changed, used in the basis of Jeremy's method of anneal value dependent on if there was no change
    double anneal_initial; // initial value for annealing, meant to replace the previously used calculation involving ANNEAL_MIN and ANNEAL_MAX with something more simple

    double mutation_rate; // fraction of new offspring to mutate
    double double_mutation_rate; // fraction of mutations to mutate two genes instead of just one
    double triple_mutation_rate; // fraction of mutations to mutate three genes instead of just one or two
    double gamma_mutate_scale; // Used in mutate(), affects the size of change for the gamma parameter values
    double tau_mutate_scale; // Used in mutate(), affects the size of change for the tau parameter values
    double coast_mutate_scale; // Used in mutate(), affects the size of change for the coast parameter values
    double triptime_mutate_scale; // Used in mutate(), affects the size of change for the triptime parameter value
    double zeta_mutate_scale; // Used in mutate(), affects the size of change for the zeta parameter value
    double beta_mutate_scale; // Used in mutate(), affects the size of change for the beta parameter value
    double alpha_mutate_scale; // Used in mutate(), affects the size of change for the alpha parameter value

    // Constructor, sets variables to default values before calling fileRead() to set them to what is in genetic.config
    geneticConstants();

    // Sets properties to what is within the config file
    // Input: Contents of genetic.config file that is in same folder as the executable file
    // Output: Properties explicitly set in the config file are set to values following 
    void geneticFileRead();
};


#include "config.cpp"

#endif
