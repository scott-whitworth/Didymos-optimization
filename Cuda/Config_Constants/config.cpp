#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include "config.h"

// Constructor first uses geneticFileRead() to set the struct's properties from a file
geneticConstants::geneticConstants() {
    // Have time_seed initially set to time(0) and it will only be changed to other value if the file has time_seed not set to -1
    this->time_seed = time(0);
    
    geneticFileRead();
}

//http://www.cplusplus.com/forum/beginner/11304/ for refesher on reading line by line
// genetic.config must be in same folder as the compiled executable
void geneticConstants::geneticFileRead() {
    std::string line;
    std::ifstream configFile;
    configFile.open("genetic.config");
    if (configFile.is_open()) {
        // Go through line by line
        while ( std::getline(configFile, line ) ) {
            // Go through if statements to find what constant variable this line refers to (look at substring prior to "=") and attempt to assign the value (after "=") to the variable
            std::string variableName = line.substr(0, line.find("=")   );
            std::string variableValue = line.substr( line.find("=") + 1 );
            
            // Assign variableValue to the appropriate variable based on variableName, with proper conversion to the right data type
            // Still need to add a check to ensure that variableValue is converted properly (if not valid, don't assign property to variableValue and note that in the terminal)
            if (variableName == "pos_threshold") {
                this->pos_threshold = std::stod(variableValue);
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
                this->triple_mutation_rate = std::stod(variableValue);
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
            else if (variableName == "time_seed" && variableValue != "-1") { // If the conifguration sets time_seed to -1 then it is left as the initially set value (time(0)) 
                this->time_seed = std::stod(variableValue);
            }
            else if (variableName == "time_seed" && variableValue == "-1") {
                std::cout << "time_seed set to time(0)\n";
            }
            else {
                std::cout << "Uknown variable '" << variableName <<"' in genetic.config!\n";
            }
        }
    }
    else {
        std::cout << "Unable to open config file!\n";
    }
}

