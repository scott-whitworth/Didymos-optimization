#include <iostream> // For cout
#include <fstream> // For file reading
#include <string>
#include <time.h> // for time(0)
#include "config.h"
#include <math.h>
#include "../constants.h" // for AU

// Constructors uses geneticFileRead() to set the struct's properties from a default config file located in same folder as executable
cudaConstants::cudaConstants() {
    FileRead("genetic.config");
    this->wet_mass = this->dry_mass + this->fuel_mass;
}
// Operates same as default, however uses configFile as address for where the config file to be used is located
cudaConstants::cudaConstants(std::string configFile) {
    FileRead(configFile);
    // Now that dry_mass and fuel_mass have been acquired, set wet_mass
    this->wet_mass = this->dry_mass + this->fuel_mass;
}


//http://www.cplusplus.com/forum/beginner/11304/ for refesher on reading line by line
void cudaConstants::FileRead(std::string fileName) {
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
                else if (variableName == "v_impact") {
                    this->v_impact = std::stod(variableValue);
                }
                else if (variableName == "rk_tol") {
                    this->rk_tol = std::stod(variableValue);
                }
                else if (variableName == "f_min") {
                    this->f_min = std::stod(variableValue);
                }
                else if (variableName == "max_numsteps") {
                    this->max_numsteps = std::stod(variableValue);
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

bool same(const cudaConstants& a, const cudaConstants& b) {
    if (a.time_seed != b.time_seed) {
        std::cout << "time_seed not equal!\n";
        return false;
    }
    else if (a.random_start != b.random_start) {
        std::cout << "random_start not equal!\n";
        return false;
    }
    else if (a.initial_start_file_address != b.initial_start_file_address) {
        std::cout << "initial file address not equal!\n";
        return false;
    }
    else if (a.pos_threshold != b.pos_threshold) {
        std::cout << "pos_threshold not equal!\n";
        return false;
    }
    else if (a.speed_threshold != b.speed_threshold) {
        std::cout << "speed_threshold not equal!\n";
        return false;
    }
    else if (a.coast_threshold != b.coast_threshold) {
        std::cout << "coast_threshold not equal!\n";
        return false;
    }
    else if (a.rk_tol != b.rk_tol) {
        std::cout << "rk_tol not equal!\n";
        return false;
    }
    else if (a.f_min != b.f_min) {
        std::cout << "f_min not equal!\n";
        return false;
    }
    else if (a.max_numsteps != b.max_numsteps) {
        std::cout << "max_numsteps not equal!\n";
        return false;
    }
    else if (a.anneal_factor != b.anneal_factor) {
        std::cout << "anneal_factor not equal!\n";
        return false;
    }
    else if (a.anneal_initial != b.anneal_initial) {
        std::cout << "anneal_initial not equal!\n";
        return false;
    }
    else if (a.change_check != b.change_check) {
        std::cout << "change_check not equal!\n";
        return false;
    }
    else if (a.v_impact != b.v_impact) {
        std::cout << "v_impact not equal!\n";
        return false;
    }
    else if (a.write_freq != b.write_freq) {
        std::cout << "write_freq not equal!\n";
        return false;
    }
    else if (a.mutation_rate != b.mutation_rate) {
        std::cout << "mutation_rate not equal!\n";
        return false;
    }
    else if (a.double_mutation_rate != b.double_mutation_rate) {
        std::cout << "double_m_rate not equal!\n";
        return false;
    }
    else if (a.triple_mutation_rate != b.triple_mutation_rate) {
        std::cout << "triple_m_rate not equal!\n";
        return false;
    }
    else if (a.gamma_mutate_scale != b.gamma_mutate_scale) {
        std::cout << "gamma_m_scale not equal!\n";
        return false;
    }
    else if (a.tau_mutate_scale != b.tau_mutate_scale) {
        std::cout << "tau_m_scale not equal!\n";
        return false;
    }
    else if (a.zeta_mutate_scale != b.zeta_mutate_scale) {
        std::cout << "zeta_m_scale not equal!\n";
        return false;
    }
    else if (a.coast_mutate_scale != b.coast_mutate_scale) {
        std::cout << "coast_m_scale not equal!\n";
        return false;
    }
    else if (a.triptime_mutate_scale != b.triptime_mutate_scale) {
        std::cout << "triptime_m_scale not equal!\n";
        return false;
    }
    else if (a.alpha_mutate_scale != b.alpha_mutate_scale) {
        std::cout << "alpha_m_scale not equal!\n";
        return false;
    }
    else if (a.beta_mutate_scale != b.beta_mutate_scale) {
        std::cout << "beta_m_scale not equal!\n";
        return false;
    }
    else if (a.thruster_type != b.thruster_type) {
        std::cout << "thruster_type not equal!\n";
        return false;
    }
    else if (a.dry_mass != b.dry_mass) {
        std::cout << "dry_mass not equal!\n";
        return false;
    }
    else if (a.fuel_mass != b.fuel_mass) {
        std::cout << "fuel_mass not equal!\n";
        return false;
    }
    else if (a.wet_mass != b.wet_mass) {
        std::cout << "wet_mass not equal!\n";
        return false;
    }
    else if (a.c3energy != b.c3energy) {
        std::cout << "c3energy not equal!\n";
        return false;
    }
    else if (a.v_escape != b.v_escape) {
        std::cout << "v_escape not equal!\n";
        return false;
    }
    else if (a.r_fin_ast != b.r_fin_ast) {
        std::cout << "r_fin_ast not equal!\n";
        return false;
    }
    else if (a.theta_fin_ast != b.theta_fin_ast) {
        std::cout << "theta_fin_ast not equal!\n";
        return false;
    }
    else if (a.z_fin_ast != b.z_fin_ast) {
        std::cout << "z_fin_ast not equal!\n";
        return false;
    }
    else if (a.vr_fin_ast != b.vr_fin_ast) {
        std::cout << "vr_fin_ast not equal!\n";
        return false;
    }
    else if (a.vtheta_fin_ast != b.vtheta_fin_ast) {
        std::cout << "vtheta_fin_ast not equal!\n";
        return false;
    }
    else if (a.vz_fin_ast != b.vz_fin_ast) {
        std::cout << "vz_fin_ast not equal!\n";
        return false;
    }
    else if (a.r_fin_earth != b.r_fin_earth) {
        std::cout << "r_fin_earth not equal!\n";
        return false;
    }
    else if (a.theta_fin_earth != b.theta_fin_earth) {
        std::cout << "theta_fin_earth not equal!\n";
        return false;
    }
    else if (a.z_fin_earth != b.z_fin_earth) {
        std::cout << "z_fin_earth not equal!\n";
        return false;
    }
    else if (a.vr_fin_earth != b.vr_fin_earth) {
        std::cout << "vr_fin_earth not equal!\n";
        return false;
    }
    else if (a.vtheta_fin_earth != b.vtheta_fin_earth) {
        std::cout << "vtheta_fin_earth not equal!\n";
        return false;
    }
    else if (a.vz_fin_earth != b.vz_fin_earth) {
        std::cout << "vz_fin_earth not equal!\n";
        return false;
    }
    else {
        std::cout << "\nThe two cudaConstants are equivalent!\n";
        return true;
    }
}


std::ostream& operator<<(std::ostream& os, const cudaConstants& object) {
    os << std::setprecision(12);
    os << "================CONFIG=DATA==============================================================\n";
    os << "time_seed: "     << object.time_seed     << "\trandom_start: "   << object.random_start      << "\t\t\taddress: "     << object.initial_start_file_address << "\n";
    os << "pos_threshold: " << object.pos_threshold << "\tspeed_threshold: "<< object.speed_threshold   << "\tcoast_threshold: " << object.coast_threshold            << "\n\n";
    
    os << "rk_tol: "        << object.rk_tol        << "\t\tf_min: "          << object.f_min           << "\t\tmax_numsteps: "  << object.max_numsteps               << "\t v_impact: " << object.v_impact << "\n";
    os << "anneal_factor: " << object.anneal_factor << "\tanneal_initial: "   << object.anneal_initial  << "\tchange_check: "    << object.change_check               << "\n";
    os << "write_freq: "    << object.write_freq    << "\t\tdisp_freq: "      << object.disp_freq       << "\t\tbest_count: "    << object.best_count                 << "\n\n";
    
    os << "mutation_rate: " << object.mutation_rate            << "\tdouble_m_rate: " << object.double_mutation_rate << "\ttriple_m_rate: "   << object.triple_mutation_rate << "\n";
    os << "gamma_m_scale: " << object.gamma_mutate_scale       << "\ttau_m_scale: "   << object.tau_mutate_scale     << "\t\tcoast_m_scale: " << object.coast_mutate_scale   << "\n";
    os << "triptime_m_scale: " << object.triptime_mutate_scale << "\talpha_m_scale: " << object.alpha_mutate_scale   << "\tbeta_m_scale: "    << object.beta_mutate_scale    << "\tzeta_m_scale: " << object.zeta_mutate_scale << "\n\n";

    os << "thruster_type: " << object.thruster_type << "\tdry_mass: " << object.dry_mass << "\t\tfuel_mass: " << object.fuel_mass << "\twet_mass: " << object.wet_mass << "\n\n";
    os << "c3energy: "      << object.c3energy      << "\tv_escape: " << object.v_escape << "\n\n";

    os << "Asteriod Info:\n";
    os << "\t R: " << object.r_fin_ast << "\t 0: " << object.theta_fin_ast << "\t Z: " << object.z_fin_ast << "\n";
    os << "\tvR: " << object.vr_fin_ast << "\tv0: " << object.vtheta_fin_ast << "\tvZ: " << object.vz_fin_ast << "\n\n";

    os << "Earth Info:\n";
    os << "\t R: " << object.r_fin_earth << "\t 0: " << object.theta_fin_earth << "\t Z: " << object.z_fin_earth << "\n";
    os << "\tvR: " << object.vr_fin_earth << "\tv0: " << object.vtheta_fin_earth << "\tvZ: " << object.vz_fin_earth << "\n\n";
    os << "===============================================================================\n";
    
    return os;
}
