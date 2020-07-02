#ifndef OUTPUT_H
#define OUTPUT_H

#include <fstream>
#include "Config_Constants/config.h"
#include "Genetic_Algorithm/individuals.h"

void trajectoryPrint(double x[], double & lastStep, int threadRank, elements<double> & yOut, thruster<double> thrust, const cudaConstants* cConstants);

void writeTrajectoryToFile(double *start, int threadRank, thruster<double> thrust, const cudaConstants* cConstants);

void progressiveAnalysis(std::ofstream & output, int rank, Individual & ind, const cudaConstants* config);

void writeIndividualToFiles(std::ofstream& ExcelOutput, std::ofstream& BinOutput, double &currentGeneration, Individual &individual, double& annealing);

void writeThrustToFiles(std::ofstream& ExcelOutput, std::ofstream& BinOutput, double &currentGeneration, Individual &individual, const cudaConstants * cConstants);

#include "output.cpp"

#endif