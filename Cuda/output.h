#ifndef OUTPUT_H
#define OUTPUT_H

#include <fstream>
#include "Config_Constants/config.h"
#include "Genetic_Algorithm/individuals.h"

// Check for error in conservation of energy due to thrust calculations
// input: 
//      time: vector of times
//      yp: vector of RK solutions (pos and vel vectors)
//      gamma: vector of in-plane thrust angles
//      tau: vector of out-of-plane thrust angles
//      lastStep: number of steps taken previously in RK method (and the length of each vector)
//      accel: vector of magnitudes of acceleration due to thrust
//      fuelSpent: vector of aggregated amount of fuel spent
//      wetMass: initial mass of the spacecraft and fuel
//      config: ptr to cudaConstants struct holding configuration values
// output:
//      errorCheck-(time_seed).bin: binary file storing vectors of time, work done, changes in mechanical energy, and avg energy between time steps
//
void errorCheck(double *time, elements<double> *yp,  double *gamma,  double *tau, int & lastStep, double *accel, double *fuelSpent, const double & wetMass, const cudaConstants* config);

// Print the trajectory (points from RK solution) to a binary file
// input: 
//      x: vector of optimized variables (initial parameters)
//      lastStep: number of steps to be taken in RK method
//      yOut: RK solution (yp) for last time step
//      thrust: thruster struct passed into rk4sys (may not be necessary since config implementation)
//      config: ptr to cudaConstants struct holding configuration values
// output:
//      orbitalMotion-(time_seed).bin: binary file storing vectors of time, yp, gamma and tau angles, thrust acceleration magnitudes, and aggregate fuel spent
//
void trajectoryPrint(double x[], double & lastStep, int threadRank, elements<double> & yOut, thruster<double> thrust, const cudaConstants* cConstants);

// Write the config parameters and optimized variables for the trajectory to a binary file
// input: 
//      start: vector of optimized variables (initial parameters)
//      threadRank: the position of this solution within the sorted pool for the current generation, offset by 1 (currently not in use; best solutions in the same run are more or less the same)
//      thrust: thruster struct passed into rk4sys (may not be necessary since config implementation)
//      config: ptr to cudaConstants struct holding configuration values
// output:
//      finalOptimization-(time_seed).bin: binary file storing vector of optimized variables, Fourier array sizes, and config values
//
void writeTrajectoryToFile(double *start, int threadRank, thruster<double> thrust, const cudaConstants* cConstants);

// Append information about the final best individual (solution) to a CSV file
// input: 
//      output: CSV file that is already open to append
//      rank: the position of this solution within the sorted pool for the current generation, offset by 1 (currently not necessary; best solutions in the same run are more or less the same)
//      ind: the individual that is being recorded
//      config: ptr to cudaConstants struct holding configuration values
// output:
//      appends data from ind to output (CSV)
//
void progressiveAnalysis(std::ofstream & output, int rank, Individual & ind, const cudaConstants* config);

// CURRENTLY NOT IN USE:
// Output functions used to analyze changes between generations within a single run

void writeIndividualToFiles(std::ofstream& ExcelOutput, std::ofstream& BinOutput, double &currentGeneration, Individual &individual, double& annealing);

// For thrust_type > 0 (only options currently are 0 and 1)
void writeThrustToFiles(std::ofstream& ExcelOutput, std::ofstream& BinOutput, double &currentGeneration, Individual &individual, const cudaConstants * cConstants);

#include "output.cpp"

#endif