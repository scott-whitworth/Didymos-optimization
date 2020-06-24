#include <iostream> // For cout
#include <fstream> // For file reading
#include <string>
#include <time.h> // for time(0)
#include "config.h"
#include <math.h>
#include "../constants.h" // for AU

// Constructors uses geneticFileRead() to set the struct's properties from a default config file located in same folder as executable
cudaConstants::cudaConstants() {
    geneticFileRead("genetic.config");
    this->wet_mass = this->dry_mass + this->fuel_mass;
}
// Operates same as default, however uses configFile as address for where the config file to be used is located
cudaConstants::cudaConstants(std::string configFile) {
    geneticFileRead(configFile);
    // Now that dry_mass and fuel_mass have been acquired, set wet_mass
    this->wet_mass = this->dry_mass + this->fuel_mass;
}


//http://www.cplusplus.com/forum/beginner/11304/ for refesher on reading line by line
void cudaConstants::geneticFileRead(std::string fileName) {
    // Use string line to hold a line read from the config file in variable configFile
    std::string line;
    std::ifstream configFile;
    configFile.open(fileName);

    if (configFile.is_open()) {
        // Go through line by line
        while ( std::getline(configFile, line ) ) {
            // If line is not empty and the line is not a comment, then the line should be a variable constant being assigned 
            if (line != "" && ( line.find("//") == std::string::npos ) && (line.find_first_of(" ") == std::string::npos) ) {
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
                else if (variableName == "coast_threshold") {
                    this->coast_threshold = std::stod(variableValue);
                }
                else if (variableName == "thruster_type") {
                    this->thruster_type = std::stoi(variableValue);
                }
                else if (variableName == "random_start") {
                    if (variableValue == "false") {
                        std::cout << "Initial generation based on file\n";
                        this->random_start = false;
                    }
                    else {
                        // If not set to false, then it is assumed the value is for true
                        std::cout << "Initial generation based on random values within set range\n";
                        this->random_start = true;
                    }
                }
                else if (variableName == "initial_start_file_address") {
                    // Assumption that the address does not need to be converted/checked
                    this->initial_start_file_address = variableValue;
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
                else if (variableName == "r_fin_ast") {
                    this->r_fin_ast = std::stod(variableValue);
                }
                else if (variableName == "theta_fin_ast") {
                    this->theta_fin_ast = std::stod(variableValue);
                }
                else if (variableName == "z_fin_ast") {
                    this->z_fin_ast = std::stod(variableValue);
                }
                else if (variableName == "vr_fin_ast") {
                    this->vr_fin_ast = std::stod(variableValue);
                }
                else if (variableName == "vtheta_fin_ast") {
                    this->vtheta_fin_ast = std::stod(variableValue);
                }
                else if (variableName == "vz_fin_ast") {
                    this->vz_fin_ast = std::stod(variableValue);
                }
                else if (variableName == "r_fin_earth") {
                    this->r_fin_earth = std::stod(variableValue);
                }
                else if (variableName == "theta_fin_earth") {
                    this->theta_fin_earth = std::stod(variableValue);
                }
                else if (variableName == "z_fin_earth") {
                    this->z_fin_earth = std::stod(variableValue);
                }
                else if (variableName == "vr_fin_earth") {
                    this->vr_fin_earth = std::stod(variableValue);
                }
                else if (variableName == "vtheta_fin_earth") {
                    this->vtheta_fin_earth = std::stod(variableValue);
                }
                else if (variableName == "vz_fin_earth") {
                    this->vz_fin_earth = std::stod(variableValue);
                }
                else if (variableName == "dry_mass") {
                   this->dry_mass = std::stod(variableValue);
                }
                else if (variableName == "fuel_mass") {
                    this->fuel_mass = std::stod(variableValue);
                }
                else if (variableName == "time_seed") { // If the conifguration sets time_seed to NONE then time_seed is set to time(0) 
                    if (variableValue != "NONE") {
                        // If variableValue is not NONE, assumption is that it is a valid double value that can be converted and used
                        this->time_seed = std::stod(variableValue);
                    }
                    else {
                        this->time_seed = time(0);
                        std::cout << "time_seed value set to time(0)\n";
                    }
                }
                else {
                    // If none of the if cases were matches, then this is some unknown variable in the config file and output this to the terminal
                    std::cout << "Unknown variable '" << variableName <<"' in " << fileName <<"!\n";
                }
            }
        }
    }
    else {
        std::cout << "Unable to open " << fileName << " file!\n";
    }
}