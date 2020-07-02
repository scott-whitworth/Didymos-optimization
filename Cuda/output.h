#ifndef OUTPUT_H
#define OUTPUT_H

#include <fstream>
#include "Config_Constants/config.h"
#include "Genetic_Algorithm/individuals.h"

void trajectoryPrint(double x[], double & lastStep, int threadRank, elements<double> & yOut, thruster<double> thrust, const cudaConstants* cConstants);

void writeTrajectoryToFile(double *start, int threadRank, thruster<double> thrust, const cudaConstants* cConstants);

void writeConfigToFile(const cudaConstants* cConstants);

void writeIndividualToFiles(std::ofstream& ExcelOutput, std::ofstream& BinOutput, double &currentGeneration, Individual &individual, double& annealing);

void writeThrustToFiles(std::ofstream& ExcelOutput, std::ofstream& BinOutput, double &currentGeneration, Individual &individual, const cudaConstants * cConstants);

#endif