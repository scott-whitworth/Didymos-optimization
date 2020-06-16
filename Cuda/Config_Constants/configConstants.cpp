#include "configReader.h"
#include <iostream>
#include <fstream>
#include <string>

// Constructor
configConstants::geneticConstants() {
    // Set this to have initial default values
    this->pos_threshold = 1.0e-10;
    this->anneal_factor = 0.25;
    this->write_freq = 100;
    this->disp_freq = 200;
    this->best_count = 10;
    this->change_check = 100;
    this->anneal_initial = 0.001;
    this->best_count = 1;

    fileRead();
}

//http://www.cplusplus.com/forum/beginner/11304/
//Currently no validation checks for converting strings to appriopriate numerica data types
void configConstants::fileRead() {
    std::string line;
    std::ifstream configFile;
    configFile.open("setup.config");
    if (configFile.is_open()) {
        while ( std::getline(configFile, line ) ) {
            std::string variableName = line.substr(0, line.find("=")   );
            std::string variableValue = line.substr( line.find("=") + 1 );

            // Assign variableValue to the appropriate variable based on variableName, with proper conversion to the right data type
            // Still need to add a check to ensure that variableValue is converted properly (if not valid, don't assign property to variableValue and note that in the terminal)
            if (variableName == "pos_threshold") {
                this->pos_threshold = std::stod(variableValue);
            }
            if (variableName == "anneal_factor") {
                this->anneal_factor = std::stod(variableValue);
            }
            if (variableName == "write_freq") {
                this->write_freq = std::stoi(variableValue);
            }
            if (variableName == "disp_freq") {
                this->disp_freq = std::stoi(variableValue);
            }
            if (variableName == "best_count") {
                this->best_count = std::stoi(variableValue);
            }
            if (variableName == "change_check") {
                this->change_check = std::stoi(variableValue);
            }
            if (variableName == "anneal_initial") {
                this->anneal_initial = std::stod(variableValue);
            }
            if (variableName = "best_count") {
                this->best_count = std::stoi(variableValue);
            }
        }
    }
    else {
        std::cout << "Unable to open config file!\n";
    }
}