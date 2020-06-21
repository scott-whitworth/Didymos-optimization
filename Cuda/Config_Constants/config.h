#ifndef GENETICCONFIG_h
#define GENETICCONFIG_h
#include <iostream>

// Structure that holds constant values related/used for the genetic algorithm that can be configured within genetic.config file
struct cudaConstants {
    double time_seed; // Seed used for randomization within optimize function, if it's set to NONE the seed is set to time(0) for genuine randomness
    double pos_threshold; // threshold for how close the spacecraft must be to the asteriod at end of its trajectory
    double speed_threshold; // threshold for desired speed the spacecraft must make
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

    // Used in thruster construction and corresponding calculations
    int thruster_type;

    double coast_threshold;
    double c3energy; // This value determines vEscape
    double v_escape; // This variable is not directly accessible in the config file as it is dependent on c3energy

    // Currently (June 19th) the following properties are not implemented
    int dry_mass;
    int wet_mass; // Wet mass is the total mass of the spacecraft (dry mass plus fuel), thinking make this derived from fuel_mass that would be in the config

    // Currently running into critical issues with this, thinking gpu usage related
    double r_fin_ast;
    double theta_fin_ast;
    double z_fin_ast;
    double vr_fin_ast;
    double vtheta_fin_ast;
    double vz_fin_ast;
    // This works fine
    double r_fin_earth;
    double theta_fin_earth;
    double z_fin_earth;
    double vr_fin_earth;
    double vtheta_fin_earth;
    double vz_fin_earth;


    // Default constructor, sets the config file path to be "genetic.config" for geneticFileRead()
    cudaConstants();

    // Constructor, accepts a string argument for the config file path
    cudaConstants(std::string configFile);

    // Sets properties to what is within the config file
    // Input: Contents of genetic.config file that is in same folder as the executable file
    // Output: Properties explicitly set in the config file are set to values following 
    void geneticFileRead(std::string fileName);
};

std::ostream& operator<<(std::ostream& os, cudaConstants& constants);

#include "config.cpp"

#endif
