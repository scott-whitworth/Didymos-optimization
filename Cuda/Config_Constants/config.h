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
    
    // Used in mutate(), affects the size of change for the respective paramater values (old code had hard-coded values)
    double gamma_mutate_scale; 
    double tau_mutate_scale; 
    double coast_mutate_scale;
    double triptime_mutate_scale;
    double zeta_mutate_scale;
    double beta_mutate_scale;
    double alpha_mutate_scale;

    // Default constructor, sets the config file path to be "genetic.config" for geneticFileRead()
    geneticConstants();

    // Constructor, accepts a string argument for the config file path
    geneticConstants(std::string configFile);

    // Sets properties to what is within the config file
    // Input: Contents of genetic.config file that is in same folder as the executable file
    // Output: Properties explicitly set in the config file are set to values following 
    void geneticFileRead(std::string fileName);
};


#include "config.cpp"

#endif
