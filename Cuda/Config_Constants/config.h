#ifndef CUDACONSTANTS_h
#define CUDACONSTANTS_h

#include <iostream>
#include <random>

// Structure that holds constant values related/used for the genetic algorithm that can be configured within genetic.config file
struct cudaConstants {
    double time_seed; // Seed used for randomization within optimize function, if it's set to NONE the seed is set to time(0) for genuine randomness
    int max_generations; // Maximum number of generations to evaluate in the genetic algortihmif not reaching a solution
    int run_count; // How many runs of the optimization algorithm want to perform with different seeds
    bool random_start; // If set to false, the initial generation has individuals initialized from a file instead of randomly generated
    bool record_mode; // If set to true, functions that record information onto files such as BestInGenerations.csv and MutateMask 
    std::string initial_start_file_address; // If random_start is false, use file_address to find what file is being used for the initial start
    
    double pos_threshold; // threshold for how close the spacecraft must be to the asteriod at end of its trajectory (AU)
    double anneal_factor; // factor by which annealing is changed when there is no change in the best individual over 100 generations

    int write_freq; // Determine how many generations between writing the progress of the best individual onto a .csv and .bin file
    int disp_freq;  // Determine how many generations between outputing contents onto the terminal screen

    int best_count; // Number of individuals that needs to be within the acceptable condition before ending the algorithm, also how many of the top individuals are recorded
    int change_check; // How often it checks for if the best individual has changed, used in the basis of Jeremy's method of anneal value dependent on if there was no change
    double anneal_initial; // initial value for annealing, meant to replace the previously used calculation involving ANNEAL_MIN and ANNEAL_MAX with something more simple

    double mutation_rate; // fraction of new offspring to mutate
    
    // Used in mutate(), affects the size of change for the respective paramater values along with annealing
    double gamma_mutate_scale; 
    double tau_mutate_scale; 
    double coast_mutate_scale;
    double triptime_mutate_scale;
    double zeta_mutate_scale;
    double beta_mutate_scale;
    double alpha_mutate_scale;

    // Used when random_start==true, is the max magnitude of the random parameter initial values (ranges from -value to +value)
    double gamma_random_start_range;
    double tau_random_start_range;
    double coast_random_start_range;
    double triptime_max; // Explicit bounds for valid triptime
    double triptime_min; // Explicit bounds for valid triptime
    double alpha_random_start_range;
    double beta_random_start_range; // For beta, only positive side of the range is used (0 to the value assigned)
    double zeta_random_start_range;

    // Used in thruster construction and corresponding calculations
    int thruster_type; // 0 is for no thruster, 1 is for NEXT ion thruster
    double dry_mass;   // Mass of the spacecraft with no fuel (kg)
    double fuel_mass;  // The mass quantity of fuel the spacecraft starts with, used to derive wet_mass
    double wet_mass;   // Wet mass is the total mass of the spacecraft (dry mass plus fuel), thinking make this derived from fuel_mass that would be in the config (kg)

    double coast_threshold; // 0 results in the thruster never coasting, 1 results in always coasting
    double c3energy;        // specific energy of spacecraft at earth escape (m^2/s^2), determines vEscape
    double v_escape;        // magnitude of velocity at earth escape  (au/s), this variable is not directly accessible in the config file as it is derived from c3energy

    // The final position and velocity of the asteriod/target at impact date
    double r_fin_ast;      // AU
    double theta_fin_ast;  // Radians
    double z_fin_ast;      // AU
    double vr_fin_ast;     // AU/s
    double vtheta_fin_ast; // AU/s
    double vz_fin_ast;     // AU/s

    // The final position and velocity of the earth at impact date to be used as reference point
    double r_fin_earth;      // AU
    double theta_fin_earth;  // Radians
    double z_fin_earth;      // AU
    double vr_fin_earth;     // AU/s
    double vtheta_fin_earth; // AU/s
    double vz_fin_earth;     // AU/s

    double v_impact; // AU/s, the official DART mission data

    double rk_tol;       // The relative/absolute (not sure which one it is) tolerance for the runge kutta algorithm
    int GuessMaxPossibleSteps; // Used as a large value to ensure adequete memory allocation in the arrays that record information in rk4sys() in output.cpp
    int cpu_numsteps; // Constant higher precision for final CPU RK method
    int min_numsteps; // Minimum time step size in runge kutta
    int max_numsteps; // Maximum time step size in runge kutta

    int num_individuals; // Number of individuals in the pool, each individual contains its own thread
    int survivor_count;  // Number of survivors selected, every pair of survivors creates 8 new individuals
    int thread_block_size;

    // Used in generating time range for Earth calculations (units in seconds)
    int timeRes;

    // Default constructor, sets the config file path to be "genetic.config" for geneticFileRead()
    cudaConstants();

    // Constructor, accepts a string argument for the config file path
    cudaConstants(std::string configFile);

    // Sets properties to what is within the config file
    // Input: File address that is used to open a text-based file and parses through to assign variables
    // Output: Properties explicitly set in the config file are set to values following equal sign, ignores comments or empty lines in files 
    void FileRead(std::string fileName);
};

// Output function to stream, with some formatting to help be more readible on terminal
std::ostream& operator<<(std::ostream& os, const cudaConstants& object);


#include "config.cpp"

#endif
