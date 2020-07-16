#ifndef OUTPUT_H
#define OUTPUT_H

#include <fstream>
#include "Config_Constants/config.h"
#include "Genetic_Algorithm/individuals.h"

// Utility function to display the currently best individual onto the terminal while the algorithm is still running
// input: Individual to be displayed (assumed to be the best individual of the pool) and the value for the current generation iterated
// output: onto the console termina, generation is displayed and best individual's posDiff, velDiff, and cost values
void terminalDisplay(Individual& individual, unsigned int currentGeneration);

// Records error in energy conservation due to thrust calculations
// Called if record_mode is true at end of optimize process
// input: time - array of times
//        yp - array of RK solutions (pos and vel components)
//        gamma - array of in-plane thrust angles
//        tau - array of out-of-plane thrust angles
//        lastStep - number of steps taken in RK
//        accel - array of magnitudes of acceleration due to thrust
//        fuelSpent - array of aggregate amounts of fuel spent
//        wetMass - initial mass of spacecraft and fuel
//        config - cudaConstants object for accessing genetic.config information
// output: errorCheck-(seed).bin recording arrays of time, work, changes in total mechanical energy, and average energy between time steps
void errorCheck(double *time, elements<double> *yp,  double *gamma,  double *tau, int & lastStep, double *accel, double *fuelSpent, const double & wetMass, const cudaConstants* config);

// Output final results of genetic algorithm
// input: x[] - array that holds parameter values of OPTIM_VARS length
//        lastStep - Used to store number of total steps as output
//        threadRank - Rank of the individual element being recorded, currently (July )
//        yOut - Used to store element informatio as output
//        thrust - Passed into rk4sys
//        cConstants - Access constants info such as target element, earth element, derive spaceCraft element, also other values such as rk_tol
// output: yOut contains final eleement information of the spacecraft
//         lastStep contains value last value for number of steps taken
void trajectoryPrint( double x[], double & lastStep, int generation, elements<double> & yOut, thruster<double> thrust, const cudaConstants* cConstants);

// Output trajectory information to finalOptimization-[time_seed].bin
// input: start - passed to trajectoryPrint, information also outputted to finalOptimization[-time_seed].bin file
//        threadRank - passed to trajectoryPrint, also 
//        thrust - passed to trajectoryPrint
//        cConstants - access time_seed for deriving file name
// output: file finalOptimization[-time_seed].bin is created that holds earth/ast/ and trajectory parameter values
void writeTrajectoryToFile(double *start, int generation, thruster<double> thrust, const cudaConstants* cConstants);

// Initialize the .csv files
// input: cConstants - to access time_seed for deriving file name conventions and also thruster type
// output: files BestInGenerations-[time_seed].csv, WorstInGenerations-[time_seed].csv, if thruster type is not NO_THRUST also BestThrustGens-[time_seed].csv & WorstThrustGens-[time_seed].csv, are given initial header row info
void initializeRecord(const cudaConstants * cConstants);

// Record progress of individual
// Called if record_mode is true at end of optimize process
// input: output - the output file stream being used
//        rank - the positional performance of the individual
//        ind - the individual object being recorded
//        config - cudaConstants object for accessing thruster_type information
// output: output file is appended information on rank, individual values/parameter information
void progressiveAnalysis(int generation, int numStep, double *start, elements<double> & yp, const cudaConstants *config);

// Utility function to observe the trend of best individual in the algorithm through the generations
// Input: Two ofstreams (one to .csv file and another to binary), current generation number, best individual, and annealing value derived to be used in next generation crossover/mutation
// Output: The two streams are appended the individual's information and anneal value
void writeIndividualToFiles(std::ofstream& ExcelOutput, std::ofstream& BinOutput, double &currentGeneration, Individual &individual, double& annealing);

// Utility function to observe the trend of best individual's thruster information in the algorithm through the generations
// Input: Two ofstreams (one to .csv file and another to binary), current generation number, best individual, and annealing value derived to be used in next generation crossover/mutation
// Output: The two streams are appended the individual's thruster information and anneal value
void writeThrustToFiles(std::ofstream& ExcelOutput, std::ofstream& BinOutput, double &currentGeneration, Individual &individual, const cudaConstants * cConstants);

// Initialize the .csv files
// input: cConstants - to access time_seed for deriving file name conventions and also thruster type
// output: files BestInGenerations-[time_seed].csv, WorstInGenerations-[time_seed].csv, if thruster type is not NO_THRUST also BestThrustGens-[time_seed].csv & WorstThrustGens-[time_seed].csv, are given initial header row info
void initializeRecord(const cudaConstants * cConstants);

// Take in the current state of the generation and appends to files
// input: cConstants - access time_seed to derive file name
//        pool - passes pool[0] to writeIndividualToFiles & writeThrustToFiles
//        generation - passed to writeIndividualToFiles & writeThrustToFiles
//        new_anneal - passed into writeIndividualToFiles
//        poolSize - used to access worst individual in pool by using index value poolSize-1
//        thrust - used in conditional statement of thrust type
// output: files BestInGenerations/WorstInGenerations have appended information using writeIndividualToFiles method
//         files BestThrustGens/WorstThurstGens have appended information using writeThrustFiles method
void recordGenerationPerformance(const cudaConstants * cConstants, Individual * pool, double generation, double new_anneal, int poolSize, thruster<double>& thrust);

// Method for doing recording information at the end of the optimization process
// input: cConstants - access record_mode, if record_mode == true then call progressiveRecord method, also passed into writeTrajectoryToFile method as well as progressiveRecord
//        pool - To access the best individual (pool[0])
//        generation - to record the generation value 
//        thrust - passed into progressiveRecord and writeTrajectoryToFile
// output: writeTrajectoryToFile is called, if in record_mode then progressiveRecord is called as well
void finalRecord(const cudaConstants* cConstants, Individual * pool, int generation, thruster<double>& thrust);

// input: cConstants - access time_seed to derive file name
// output: mutateFile[time_seed].csv is given a header row, now ready to be used for progressiveRecord()
void setMutateFile(const cudaConstants* cConstants);

// method that stores information of launchCon of timeRes*24 resolution
// input: cConstants - access time range and resolution info on launchCon
//        launchCon - access elements of earth 
// output: EarthCheckValues.csv is created and holds rows of element info on earth with timeStamp on each row
void recordEarthData(const cudaConstants * cConstants);

// Takes in a pool and records the parameter info on all individuals
// input: cConstants - to access time_seed in deriving file name
//        pool - holds all the individuals to be stored
//        poolSize - to use in iterating through the pool
//        generation - used in deriving file name
// output: file generation#[generation]-[time_seed].csv is created with each row holding parameter values of individuals
void recordAllIndividuals(const cudaConstants * cConstants, Individual * pool, int poolSize, int generation);

#include "output.cpp"

#endif