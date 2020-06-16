#ifndef CONFIGREADER_h
#define CONFIGREADER_h


// configConstants is a struct to hold the various constant values used in the program that we would want to be able to change without needing to recompile everything
// OLD APPROACH -> Currently working on utilizing enumerations as way of storing values
struct configConstants {
    // Constants we want to change that are currently found in gaConstants.h
    double pos_threshold;
    double anneal_initial;
    double anneal_factor;
    int write_freq;
    int disp_freq;
    int best_count;
    int change_check;

    // Constants we want to change that are currently found in constants.h
    double coast_threshold;
    double thruster_id;
    // Elements of asteriod and earth?
    // Spacecraft constants? (Andrew - thinking we should have dry_mass and fuel_mass)

    // Constructor, sets variables to default values before calling fileRead() to set them to what is in setup.config
    geneticConstants();

    // Opens setup.config and sets properties to what is within the file
    // Assumption is that setup.config is in the same folder as the executable file
    fileRead();
}

#include "configConstants.cpp"

#endif