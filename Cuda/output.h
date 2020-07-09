#ifndef OUTPUT_H
#define OUTPUT_H

#include <fstream>
#include "Config_Constants/config.h"
#include "Genetic_Algorithm/individuals.h"

double calc_work(double time, elements<double> yp, double gamma, double tau, double accel, double mass);

double calc_dE(double time, elements<double> yp, double mass);

void errorCheck(double *time, elements<double> *yp,  double *gamma,  double *tau, int & lastStep, double *accel, double *fuelSpent, const double & wetMass, const cudaConstants* config);

void trajectoryPrint(double x[], double & lastStep, int threadRank, elements<double> & yOut, thruster<double> thrust, const cudaConstants* cConstants);

void writeTrajectoryToFile(double *start, int threadRank, thruster<double> thrust, const cudaConstants* cConstants);

void progressiveAnalysis(std::ofstream & output, int rank, Individual & ind, const cudaConstants* config);

void writeIndividualToFiles(std::ofstream& ExcelOutput, std::ofstream& BinOutput, double &currentGeneration, Individual &individual, double& annealing);

void writeThrustToFiles(std::ofstream& ExcelOutput, std::ofstream& BinOutput, double &currentGeneration, Individual &individual, const cudaConstants * cConstants);

#include "output.cpp"

#endif