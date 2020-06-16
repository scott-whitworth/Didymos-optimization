#include <string>
#include <fstream>
#include <iostream>
// WIP - Attempt to have the configuration constants be enumerated constants
namespace config {
    enum class geneticConstants {
        pos_threshold,
        anneal_initial,
        anneal_factor,
        write_freq,
        disp_freq,
        best_count,
        change_check
    };

    void geneticFileRead() {
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
                    geneticConstants::pos_threshold = std::stod(variableValue);
                }
                if (variableName == "anneal_factor") {
                    anneal_factor = std::stod(variableValue);
                }
                if (variableName == "write_freq") {
                    write_freq = std::stoi(variableValue);
                }
                if (variableName == "disp_freq") {
                    disp_freq = std::stoi(variableValue);
                }
                if (variableName == "best_count") {
                    best_count = std::stoi(variableValue);
                }
                if (variableName == "change_check") {
                    change_check = std::stoi(variableValue);
                }
                if (variableName == "anneal_initial") {
                    anneal_initial = std::stod(variableValue);
                }
                if (variableName = "best_count") {
                    best_count = std::stoi(variableValue);
                }
            }
        }
        else {
            std::cout << "Unable to open config file!\n";
        }

    }
}