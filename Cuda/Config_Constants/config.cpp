#include <iostream>
#include <fstream>
#include <string>

#include "config.h"

geneticConstants::geneticConstants() {
    // Set this to have initial values in case the file doesn't have them set to anything
    this->pos_threshold = 1.0e-10;
    
    this->anneal_factor = 0.25;
    this->anneal_initial = 0.001;
    
    this->write_freq = 100;
    this->disp_freq = 200;
    this->best_count = 10;
    this->change_check = 100;
    
    this->mutation_rate = 0.15;
    this->double_mutation_rate = 0.2;
    this->triple_mutation_rate = 0.05;

    this->gamma_mutate_scale=10.0;
    this->tau_mutate_scale=10.0;
    this->coast_mutate_scale=10.0;
    this->triptime_mutate_scale=0.5;
    this->zeta_mutate_scale=1.57;
    this->beta_mutate_scale=3.14;
    this->alpha_mutate_scale=3.14;
    
    // NOT COMPLETE ENSURE ALL ARE ADDED
    geneticFileRead();
}

//http://www.cplusplus.com/forum/beginner/11304/
void geneticConstants::geneticFileRead() {
    std::string line;
    std::ifstream configFile;
    configFile.open("genetic.config");
    if (configFile.is_open()) {
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
            else if (variableName == "double_mutation_rate") {
                this->double_mutation_rate = std::stod(variableValue);
            }
            else if (variableName == "triple_mutation_rate") {
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
            else {
                std::cout << "Uknown variableName in genetic.config!\n";
            }
        }
    }
    else {
        std::cout << "Unable to open config file!\n";
    }
}

