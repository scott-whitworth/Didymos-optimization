#include "output.h"
#include "constants.h"
#include "Thrust_files/thruster.h"
#include <string>
#include <iomanip>

// Output final results of genetic algorithm
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------
void trajectoryPrint( double x[], double & lastStep, int threadRank, elements<double> & yOut, thruster<double> thrust, const cudaConstants* cConstants) {
  /*set the asteroid and inital conditions for the earth and spacecraft:
  constructor takes in radial position(au), angluar position(rad), axial position(au),
  radial velocity(au/s), tangential velocity(au/s), axial velocity(au/s)*/

  // setting landing conditions of the asteroid (Sept 30, 2022)
  elements<double> asteroid = elements<double>(cConstants->r_fin_ast, cConstants->theta_fin_ast, cConstants->z_fin_ast, cConstants->vr_fin_ast, cConstants->vtheta_fin_ast, cConstants->vz_fin_ast);

  // setting initial conditions of earth based off of the impact date (Sept 30, 2022) minus the trip time (optimized).
  elements<double> earth =  launchCon->getCondition(x[TRIPTIME_OFFSET]);
  
  // setting initial conditions of the spacecraft
  elements<double> spaceCraft = elements<double>(earth.r+ESOI*cos(x[ALPHA_OFFSET]),
                                                 earth.theta + asin(sin(M_PI-x[ALPHA_OFFSET])*ESOI/earth.r),
                                                 earth.z,
                                                 earth.vr + cos(x[ZETA_OFFSET])*sin(x[BETA_OFFSET])*cConstants->v_escape,
                                                 earth.vtheta + cos(x[ZETA_OFFSET])*cos(x[BETA_OFFSET])*cConstants->v_escape,
                                                 earth.vz + sin(x[ZETA_OFFSET])*cConstants->v_escape);

  // setting time parameters
  double timeInitial=0; 
  double timeFinal=orbitalPeriod; // Orbital period of asteroid(s)
  double deltaT; // time step
  int numSteps = 5000; // initial guess for the number of time steps, guess for the memory allocated 
  deltaT = (timeFinal-timeInitial) / cConstants->max_numsteps; // initial guess for time step, small is preferable

  // setup of thrust angle calculations based off of optimized coefficients
  coefficients<double> coeff;
  initCoefficient(x,coeff, cConstants);
  // Assigning coast threshold (now done in coefficients because is a constant)

  // Assigning wetMass
  double wetMass = cConstants->wet_mass;
//  double wetMass = WET_MASS;
  // setting Runge-Kutta tolerance
  double absTol = cConstants->rk_tol;
  // set optmization minimum
  // double Fmin = cConstants->f_min;

  // Initialize memory for the solution vector of the dependant solution
  elements<double>* yp;
  yp = new elements<double>[numSteps];
  
  double *times, *gamma, *tau, *accel_output, *fuelSpent;
  times = new double[numSteps]; // Initialize memory for time array
  gamma = new double[numSteps]; // Initialize memory for gamma array
  tau = new double[numSteps]; // Initialize memory for tau array
  accel_output = new double[numSteps]; // Initialize memory for acceleration array
  fuelSpent = new double[numSteps];  // Initialize memory for fuelSpent array

  double accel; // Initialize memory for  acceleration

  // used to track the cost function throughout a run via output and outputs to a binary
  int lastStepInt;

  rk4sys(timeInitial, x[TRIPTIME_OFFSET] , times, spaceCraft, deltaT, yp, absTol, coeff, accel, gamma, tau, lastStepInt, accel_output, fuelSpent, wetMass, thrust, cConstants);

  lastStep = lastStepInt;

  // gets the final y values of the spacecrafts for the cost function.
  yOut = yp[(int)lastStep];

  std::ofstream output;
  double seed = cConstants->time_seed;
  output.open("orbitalMotion-"+std::to_string(static_cast<int>(seed))+".bin", std::ios::binary);
  // output.open("orbitalMotion-"+std::to_string(static_cast<int>(seed))+"-"+std::to_string(threadRank)+".bin", std::ios::binary);
  for(int i = 0; i <= lastStep; i++) {
    //output << yp[i];
    output.write((char*)&yp[i], sizeof (elements<double>));
    output.write((char*)&times[i], sizeof (double));
    output.write((char*)&gamma[i], sizeof (double));
    output.write((char*)&tau[i], sizeof (double));
    output.write((char*)&accel_output[i], sizeof (double));
    output.write((char*)&fuelSpent[i], sizeof (double));
  }
  output.close();
  
  // cleaning up dynamic yp, time, gamma, and tau.
  delete [] yp;
  delete [] times;
  delete [] gamma;
  delete [] tau;
}

void writeTrajectoryToFile(double *start, int threadRank, thruster<double> thrust, const cudaConstants* cConstants) {
    elements<double> yp;
    double numStep = 0;
    trajectoryPrint(start, numStep, threadRank, yp, thrust, cConstants);

    //writes final optimization values to a seperate file
    std::ofstream output;
    double seed = cConstants->time_seed;
    output.open("finalOptimization-"+std::to_string(static_cast<int>(seed))+".bin", std::ios::binary);
    // output.open ("finalOptimization-"+std::to_string(static_cast<int>(seed))+"-"+std::to_string(threadRank)+".bin", std::ios::binary);

    output.write((char*)&cConstants->r_fin_ast, sizeof(double));
    output.write((char*)&cConstants->theta_fin_ast, sizeof(double));
    output.write((char*)&cConstants->z_fin_ast, sizeof(double));
    output.write((char*)&cConstants->vr_fin_ast, sizeof(double));
    output.write((char*)&cConstants->vtheta_fin_ast, sizeof(double));
    output.write((char*)&cConstants->vz_fin_ast, sizeof(double));
    output.write((char*)&cConstants->r_fin_earth, sizeof(double));
    output.write((char*)&cConstants->theta_fin_earth, sizeof(double));
    output.write((char*)&cConstants->z_fin_earth, sizeof(double));
    output.write((char*)&cConstants->vr_fin_earth, sizeof(double));
    output.write((char*)&cConstants->vtheta_fin_earth, sizeof(double));
    output.write((char*)&cConstants->vz_fin_earth, sizeof(double));
    output.write((char*)&cConstants->coast_threshold, sizeof(double));
    double gsize = GAMMA_ARRAY_SIZE, tsize = TAU_ARRAY_SIZE, csize = COAST_ARRAY_SIZE;
    output.write((char*)&gsize, sizeof(double));
    output.write((char*)&tsize, sizeof(double));
    output.write((char*)&csize, sizeof(double));

    for (int j = 0; j < OPTIM_VARS; j++) {
      output.write((char*)&start[j], sizeof (double));
    }
    
    output.write((char*)&numStep, sizeof (double));

  output.close();
}

void progressiveAnalysis(std::ofstream & output, int rank, Individual & ind, const cudaConstants* config) {
    output << rank << ',' << ind.posDiff << ',' << ind.velDiff << ',' << ind.startParams.tripTime << ',';
    output << ind.startParams.alpha << ',' << ind.startParams.beta << ',' << ind.startParams.zeta << ',';
    if (config->thruster_type) {
      for (int i = 0; i < GAMMA_ARRAY_SIZE; i++) {
        output << ind.startParams.coeff.gamma[i] << ',';
      }
      for (int i = 0; i < TAU_ARRAY_SIZE; i++) {
        output << ind.startParams.coeff.tau[i] << ',';
      }
      for (int i = 0; i < COAST_ARRAY_SIZE; i++) {
        output << ind.startParams.coeff.coast[i] << ',';
      }
    }
    output << std::endl;
}

// Initialize the .csv files
void initializeRecord(const cudaConstants * cConstants) {
    // setup output of generation results over time onto a .csv file
    std::ofstream bestExcel;
    bestExcel.open("BestInGenerations-"+ std::to_string(cConstants->time_seed)+".csv");
    // Set first row in the file be a header for the columns
    bestExcel << "Gen #" << "," << "posDiff" << "," << "velDiff" << "," 
                                   << "rFinal" << "," << "thetaFinal" << "," << "zFinal" << "," << "vrFinal" << "," << "vthetaFinal" << "," << "vzFinal" << ","
                                   << "rInitial" << "," << "thetaInitial" << "," << "zInitial" << ","<< "vrInitial" << "," << "vthetaInitial" << "," << "vzInitial" << ","
                                   << "alpha" << "," << "beta" << "," << "zeta" << "," << "anneal" << "," << "tripTime" << "\n";
    bestExcel.close();
    std::ofstream worstExcel;
    worstExcel.open("WorstInGenerations-"+ std::to_string(cConstants->time_seed)+".csv");
    // Set first row in the file be a header for the columns
    worstExcel << "Gen #" << "," << "posDiff" << "," << "velDiff" << "," 
                                << "rFinal" << "," << "thetaFinal" << "," << "zFinal" << "," << "vrFinal" << "," << "vthetaFinal" << "," << "vzFinal" << ","
                                << "rInitial" << "," << "thetaInitial" << "," << "zInitial" << ","<< "vrInitial" << "," << "vthetaInitial" << "," << "vzInitial" << ","
                                << "alpha" << "," << "beta" << "," << "zeta" << "," << "anneal" << "," << "tripTime" << "\n";
    worstExcel.close();

    if (cConstants->thruster_type != thruster<double>::NO_THRUST) {
        std::ofstream thrustBestExcel, thrustWorstExcel;

        thrustBestExcel.open("BestThrustGens-"+ std::to_string(cConstants->time_seed)+".csv");
        thrustBestExcel << "gen,gamma0,gamma1,gamma2,gamma3,gamma4,gamma5,gamma6,tau0,tau1,tau2,coast0,coast1,coast2,coast3,coast4";

        thrustWorstExcel.open("WorstThrustGens-"+ std::to_string(cConstants->time_seed)+".csv");
        thrustWorstExcel << "gen,gamma0,gamma1,gamma2,gamma3,gamma4,gamma5,gamma6,tau0,tau1,tau2,coast0,coast1,coast2,coast3,coast4";
    }
    
}

// Take in the current state of the generation and appends to files
void recordGenerationPerformance(const cudaConstants * cConstants, Individual * pool, double generation, double new_anneal, int poolSize) {
  std::ofstream bestExcel, bestBin;
  bestExcel.open("BestInGenerations-"+ std::to_string(cConstants->time_seed)+".csv", std::ios_base::app);
  bestBin.open("BestInGenerations-"+ std::to_string(cConstants->time_seed)+".bin", std::ios_base::app);

  writeIndividualToFiles(bestExcel, bestBin, generation, pool[0], new_anneal);

  bestExcel.close();
  bestBin.close();

  std::ofstream worstExcel, worstBin;
  worstExcel.open("WorstInGenerations-"+ std::to_string(cConstants->time_seed)+".csv", std::ios_base::app);
  worstBin.open("WorstInGenerations-"+ std::to_string(cConstants->time_seed)+".bin", std::ios_base::app);

  writeIndividualToFiles(worstExcel, worstBin, generation, pool[poolSize-1], new_anneal);

  worstExcel.close();
  worstBin.close();

  if (cConstants->thruster_type != thruster<double>::NO_THRUST) {
    std::ofstream bestThrusterExcel, bestThrusterBin;
    bestThrusterExcel.open("BestThrustGens-"+ std::to_string(cConstants->time_seed)+".csv", std::ios_base::app);
    bestThrusterBin.open("BestThurstGens-"+ std::to_string(cConstants->time_seed)+".bin", std::ios_base::app);
    
    writeThrustToFiles(bestThrusterExcel, bestThrusterBin, generation, pool[0], cConstants);
    
    bestExcel.close();
    bestBin.close();

    std::ofstream worstThrusterExcel, worstThrusterBin;
    worstExcel.open("WorstThrustGens-"+ std::to_string(cConstants->time_seed)+".csv", std::ios_base::app);
    worstBin.open("WorstThrustGens-"+ std::to_string(cConstants->time_seed)+".bin", std::ios_base::app);

    writeThrustToFiles(worstThrusterExcel, worstThrusterBin, generation, pool[poolSize-1], cConstants);
    
    worstThrusterExcel.close();
    worstThrusterBin.close();
  }

}



// Utility function to observe the trend of best individual in the algorithm through the generations
// Input: Two ofstreams (one to .csv file and another to binary), current generation number, best individual, and annealing value derived to be used in next generation crossover/mutation
// Output: The two streams are appended the individual's information and anneal value
void writeIndividualToFiles(std::ofstream& ExcelOutput, std::ofstream& BinOutput, double &currentGeneration, Individual &individual, double& annealing ) {
    // Output the information to excel spreadsheet
    ExcelOutput << currentGeneration << ','
                << individual.posDiff << ',' << individual.velDiff << ',' // The positional and velocity difference
                << individual.finalPos.r << ',' << individual.finalPos.theta << ',' << individual.finalPos.z << ',' // Final position
                << individual.finalPos.vr << ',' << individual.finalPos.vtheta << ',' << individual.finalPos.vz << ',' // Final velocity
                << individual.startParams.y0.r << ',' << individual.startParams.y0.theta << ',' << individual.startParams.y0.z << ',' // Starting position
                << individual.startParams.y0.vr << ',' << individual.startParams.y0.vtheta << ',' << individual.startParams.y0.vz << ',' // Starting velocity
                << individual.startParams.alpha << ',' << individual.startParams.beta << ',' << individual.startParams.zeta << ',' // alpha, beta, zeta
                << annealing << "," << individual.startParams.tripTime << std::endl; // Annealing value for next generation and triptime (in that order to maintain continuity with bin file)
 
    // Output the information to binary file for use in the MATLAB code, line breaks and spaces added to help with readibility
    BinOutput.write( (char*)& currentGeneration, sizeof(double));
    // posDiff and velDiff
    BinOutput.write( (char*)& individual.posDiff, sizeof(double));
    BinOutput.write( (char*)& individual.velDiff, sizeof(double));
    // Position and velocity information
    BinOutput.write( (char*)& individual.finalPos.r,            sizeof(double));
    BinOutput.write( (char*)& individual.finalPos.theta,        sizeof(double));
    BinOutput.write( (char*)& individual.finalPos.z,            sizeof(double));
    BinOutput.write( (char*)& individual.finalPos.vr,           sizeof(double));
    BinOutput.write( (char*)& individual.finalPos.vtheta,       sizeof(double));
    BinOutput.write( (char*)& individual.finalPos.vz,           sizeof(double));
    BinOutput.write( (char*)& individual.startParams.y0.r,      sizeof(double));
    BinOutput.write( (char*)& individual.startParams.y0.theta,  sizeof(double));
    BinOutput.write( (char*)& individual.startParams.y0.z,      sizeof(double));
    BinOutput.write( (char*)& individual.startParams.y0.vr,     sizeof(double));
    BinOutput.write( (char*)& individual.startParams.y0.vtheta, sizeof(double));
    BinOutput.write( (char*)& individual.startParams.y0.vz,     sizeof(double));
    // Alpha, Beta, Zeta, Annealing, Triptime
    BinOutput.write( (char*)& individual.startParams.alpha,  sizeof(double));
    BinOutput.write( (char*)& individual.startParams.beta,   sizeof(double));
    BinOutput.write( (char*)& individual.startParams.zeta,   sizeof(double));
    BinOutput.write((char*)& annealing, sizeof(double));
    BinOutput.write((char*)& individual.startParams.tripTime, sizeof(double));
}

void writeThrustToFiles(std::ofstream& ExcelOutput, std::ofstream& BinOutput, double &currentGeneration, Individual &individual, const cudaConstants * cConstants) {
    ExcelOutput << currentGeneration << ',';
    for (int i = 0; i < GAMMA_ARRAY_SIZE; i++) {
        ExcelOutput << individual.startParams.coeff.gamma[i] << ',';
    }
    for (int i = 0; i < TAU_ARRAY_SIZE; i++) {
        ExcelOutput << individual.startParams.coeff.tau[i] << ',';
    }
    for (int i = 0; i < COAST_ARRAY_SIZE; i++) {
        ExcelOutput << individual.startParams.coeff.coast[i] << ',';
    }
    ExcelOutput << std::endl;

    BinOutput.write((char*)&currentGeneration, sizeof(double));
    for (int i = 0; i < GAMMA_ARRAY_SIZE; i++) {
        BinOutput.write((char*)&individual.startParams.coeff.gamma[i], sizeof(double));
    }
    for (int i = 0; i < TAU_ARRAY_SIZE; i++) {
        BinOutput.write((char*)&individual.startParams.coeff.tau[i], sizeof(double));
    }
    for (int i = 0; i < COAST_ARRAY_SIZE; i++) {
        BinOutput.write((char*)&individual.startParams.coeff.coast[i], sizeof(double));
    }
}