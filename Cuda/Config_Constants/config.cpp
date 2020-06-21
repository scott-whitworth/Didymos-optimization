#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include "config.h"
#include <math.h>
#include "../constants.h"

using namespace std;
using std::stod;
// Constructors uses geneticFileRead() to set the struct's properties from a default config file located in same folder as executable
cudaConstants::cudaConstants() {
    geneticFileRead("genetic.config");
}
// Operates same as default, however uses configFile as address for where the config file to be used is located
cudaConstants::cudaConstants(std::string configFile) {
    geneticFileRead(configFile);
}


//http://www.cplusplus.com/forum/beginner/11304/ for refesher on reading line by line
void cudaConstants::geneticFileRead(std::string fileName) {

    std::string line;
    std::ifstream configFile;
    configFile.open(fileName);
    if (configFile.is_open()) {
        // Go through line by line
        while ( std::getline(configFile, line ) ) {
            // If line is not empty and the line is not a comment, then the line should be a variable constant being assigned 
            if (line != "" && ( line.find_first_of("//") == std::string::npos ) && (line.find_first_of(" ") == std::string::npos) ) {
                // Go through if statements to find what constant variable this line refers to (look at substring prior to "=") and attempt to assign the value (after "=") to the variable
                std::string variableName = line.substr(0, line.find("=")   );
                std::string variableValue = line.substr( line.find("=") + 1); // Currently reads variableValue to end of the line string, may want to implement an end position to allow comments appends in the file format

                // Assign variableValue to the appropriate variable based on variableName, with proper conversion to the right data type
                // Still need to add a check to ensure that variableValue is converted properly (if not valid, don't assign property to variableValue and note that in the terminal)
                if (variableName == "pos_threshold") {
                    this->pos_threshold = std::stod(variableValue);
                }
                else if (variableName == "speed_threshold") {
                    this->speed_threshold = std::stod(variableValue);
                }
                else if (variableName == "thruster_type") {
                    this->thruster_type = std::stoi(variableValue);
                }
                else if (variableName == "anneal_factor") {
                    this->anneal_factor = std::stod(variableValue);
                }
                else if (variableName == "write_freq") {
                    this->write_freq = std::stoi(variableValue);
                }
                else if (variableName == "disp_freq") {
                    this->disp_freq = std::stoi(variableValue);
                }
                else if (variableName == "best_count") {
                    this->best_count = std::stoi(variableValue);
                }
                else if (variableName == "change_check") {
                    this->change_check = std::stoi(variableValue);
                }
                else if (variableName == "anneal_initial") {
                    this->anneal_initial = std::stod(variableValue);
                }
                else if (variableName == "mutation_rate") {
                    this->mutation_rate = std::stod(variableValue);
                }
                else if (variableName == "double_mutate_rate") {
                    this->double_mutation_rate = std::stod(variableValue);
                }
                else if (variableName == "triple_mutate_rate") {
                    this->triple_mutation_rate = std::stod(variableValue);
                }
                else if (variableName == "gamma_mutate_scale") {
                    this->gamma_mutate_scale = std::stod(variableValue);
                }
                else if (variableName == "tau_mutate_scale") {
                    this->tau_mutate_scale = std::stod(variableValue);
                }
                else if (variableName == "coast_mutate_scale") {
                    this->coast_mutate_scale = std::stod(variableValue);
                }
                else if (variableName == "triptime_mutate_scale") {
                    this->triptime_mutate_scale = std::stod(variableValue);
                }
                else if (variableName == "zeta_mutate_scale") {
                    this->zeta_mutate_scale = std::stod(variableValue);
                }
                else if (variableName == "beta_mutate_scale") {
                    this->beta_mutate_scale = std::stod(variableValue);
                }
                else if (variableName == "alpha_mutate_scale") {
                    this->alpha_mutate_scale = std::stod(variableValue);
                }
                else if (variableName == "coast_threshold") {
                    this->coast_threshold = std::stod(variableValue);
                }
                else if (variableName == "c3energy") { // Determine's not just c3energy, but also v_escape
                    this->c3energy = std::stod(variableValue);
                    this->v_escape = sqrt(this->c3energy)/AU;
                }
//                else if (variableName == "r_fin_ast") {
//                    this->r_fin_ast = std::stod(variableValue);
//                }
//                else if (variableName == "theta_fin_ast") {
//                    this->theta_fin_ast = std::stod(variableValue);
//                }
//                else if (variableName == "z_fin_ast") {
//                    this->z_fin_ast = std::stod(variableValue);
//                }
//                else if (variableName == "vr_fin_ast") {
//                    this->vr_fin_ast = std::stod(variableValue);
//                }
//                else if (variableName == "vtheta_fin_ast") {
//                    this->vtheta_fin_ast = std::stod(variableValue);
//                }
//                else if (variableName == "vz_fin_ast") {
//                    this->vz_fin_ast = std::stod(variableValue);
//                }
//                else if (variableName == "r_fin_earth") {
//                    this->r_fin_earth = std::stod(variableValue);
//                }
//                else if (variableName == "theta_fin_earth") {
//                    this->theta_fin_earth = std::stod(variableValue);
//                }
//                else if (variableName == "z_fin_earth") {
//                    this->z_fin_earth = std::stod(variableValue);
//                }
//                else if (variableName == "vr_fin_earth") {
//                    this->vr_fin_earth = std::stod(variableValue);
//                }
//                else if (variableName == "vtheta_fin_earth") {
//                    this->vtheta_fin_earth = std::stod(variableValue);
//                }
//                else if (variableName == "vz_fin_earth") {
//                    this->vz_fin_earth = std::stod(variableValue);
//                }
                else if (variableName == "dry_mass") {
                    this->dry_mass = std::stoi(variableValue);
                }
                else if (variableName == "wet_mass") {
                    this->wet_mass = std::stoi(variableValue);
                }
                else if (variableName == "time_seed" && variableValue != "NONE") { // If the conifguration sets time_seed to NONE then time_seed is set to (time(0)) 
                    this->time_seed = std::stod(variableValue);
                }
                else if (variableName == "time_seed" && variableValue == "NONE") {
                    this->time_seed = time(0);
                    std::cout << "time_seed set to time(0)\n";
                }
                else {
                    std::cout << "Unknown variable '" << variableName <<"' in " << fileName <<"!\n";
                }
            }
        }
    }
    else {
        std::cout << "Unable to open " << fileName << " file!\n";
    }
}



std::ostream& operator<<(std::ostream& os, const cudaConstants& constants ) {
    os << "-Config File Reading Results-\n";
    os << "time_seed:" << constants.time_seed;
    os << "\npos_threshold:" << constants.pos_threshold;
    os << "\nspeed_threshold: " << constants.speed_threshold;
    os << "\nanneal_factor: " << constants.anneal_factor;
    os << "\nwrite_freq: " << constants.write_freq;
    os << "\ndisp_freq: " << constants.disp_freq;
    os << "\nbest_count: " << constants.best_count;
    os << "\nchange_check: " << constants.change_check;
    os << "\nanneal_initial: " << constants.anneal_initial;
    os << "\nmutation_rate: " << constants.mutation_rate; 
    os << "\ndouble_mutation_rate: " << constants.double_mutation_rate;
    os << "\ntriple_mutation_rate: " << constants.triple_mutation_rate;
    os << "\ngamma_mutate_scale: " << constants.gamma_mutate_scale; 
    os << "\ntau_mutate_scale: " << constants.tau_mutate_scale;
    os << "\ncoast_mutate_scale: " << constants.coast_mutate_scale;
    os << "\ntriptime_mutate_scale: " << constants.triptime_mutate_scale;
    os << "\nzeta_mutate_scale: " << constants.zeta_mutate_scale;
    os << "\nbeta_mutate_scale: "<< constants.beta_mutate_scale;
    os << "\nalpha_mutate_scale: " << constants.alpha_mutate_scale;
    os << "\nthruster_type: " << constants.thruster_type;
    os << "\ncoast_threshold: " << constants.coast_threshold;
    os << "\nc3energy: " << constants.c3energy;/*
    os << "\nr_fin_ast: " << constants.r_fin_ast;
    os << "\ntheta_fin_ast: " << constants.theta_fin_ast;
    os << "\nz_fin_ast: " << constants.z_fin_ast;
    os << "\nvr_fin_ast: " << constants.vr_fin_ast;
    os << "\nvtheta_fin_ast: " << constants.vtheta_fin_ast;
    os << "\nvz_fin_ast: " << constants.vz_fin_ast;
    os << "\nr_fin_earth: " << constants.r_fin_earth;
    os << "\ntheta_fin_earth: " << constants.theta_fin_earth;
    os << "\nz_fin_earth: " << constants.z_fin_earth;
    os << "\nvr_fin_earth: " << constants.vr_fin_earth;
    os << "\nvtheta_fin_earth: " << constants.vtheta_fin_earth;
    os << "\nvz_fin_earth: " << constants.vz_fin_earth;*/
    os << "\ndry_mass: " << constants.dry_mass;
    os << "\nwet_mass: " << constants.wet_mass;
    return os;
}