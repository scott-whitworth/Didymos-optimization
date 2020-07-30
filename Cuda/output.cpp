#include "output.h"
#include "constants.h"
#include "Thrust_files/thruster.h"
#include <string>
#include <iomanip>
#include "math.h"

// Utility function to display the currently best individual onto the terminal while the algorithm is still running
// input: Individual to be displayed (assumed to be the best individual of the pool) and the value for the current generation iterated
// output: onto the console termina, generation is displayed and best individual's posDiff, velDiff, and cost values
void terminalDisplay(Individual& individual, unsigned int currentGeneration) {
    std::cout << "\nGeneration: " << currentGeneration << std::endl;
    std::cout << "Best individual:" << std::endl;
    std::cout << "\tposDiff: " << individual.posDiff << std::endl;
    std::cout << "\tvelDiff: " << individual.velDiff << std::endl;
    std::cout << "\tcost: "    << individual.cost << std::endl;
}

// input: cConstants - access time_seed to derive file name
// output: mutateFile[time_seed].csv is given a header row, now ready to be used for progressiveRecord()
void setMutateFile(const cudaConstants* cConstants) { 
  std::ofstream mutateFile;
  int seed = cConstants->time_seed;
  mutateFile.open("mutateFile-" + std::to_string(seed) + ".csv", std::ios_base::app);

  mutateFile << "gen, anneal, genesToMutate,";
  for (int i = 0; i < GAMMA_ARRAY_SIZE; i++) {
    mutateFile << "gamma" << i << ",";
  }
  for (int i = 0; i < TAU_ARRAY_SIZE; i++) {
    mutateFile << "tau" << i << ",";
  }
  for (int i = 0; i < COAST_ARRAY_SIZE; i++) {
    mutateFile << "coast" << i << ",";
  }
  mutateFile << "alpha,beta,zeta,tripTime\n";

  mutateFile.close();
}

void errorCheck(double *time, elements<double> *yp,  double *gamma,  double *tau, int & lastStep, double *accel, double *fuelSpent, const double & wetMass, double *work, double *dE, double *Etot_avg, const cudaConstants* config) {
  double *mass = new double[lastStep];
  double *Etot = new double[lastStep];
  
  for (int i = 0; i < lastStep; i++) {
    mass[i] = wetMass - fuelSpent[i];
    // W = (F.v)dt (less accurate; implementation uses old index format without averages)
    // work[i] = mass[i] * accel[i] * time[i] * sqrt(pow(sin(gamma[i])*cos(tau[i])*yp[i+1].vr, 2) + pow(cos(gamma[i])*cos(tau[i])*yp[i].vtheta, 2) + pow(sin(tau[i])*yp[i].vz, 2)) / pow(AU,2);
    Etot[i] = mass[i] * ((pow(yp[i].vr,2) + pow(yp[i].vtheta,2) + pow(yp[i].vz,2))/ 2 - constG * massSun / yp[i].r) / pow(AU,2);
    if (i) {
      // W = F.dL (more accurate)
      work[i] = (mass[i] + mass[i-1])/2 * (accel[i] + accel[i-1])/2 * ((sin((gamma[i] + gamma[i-1])/2)*cos((tau[i] + tau[i-1])/2)*(yp[i].r - yp[i-1].r)) + (cos((gamma[i] + gamma[i-1])/2)*cos((tau[i] + tau[i-1])/2)*(yp[i].r + yp[i-1].r)/2*(yp[i].theta - yp[i-1].theta)) + (sin((tau[i] + tau[i-1])/2)*(yp[i].z - yp[i-1].z))) / pow(AU,2);
      dE[i] = Etot[i] - Etot[i-1];
      Etot_avg[i] = (Etot[i] + Etot[i-1])/2;
    }
  }
  work[0] = dE[0] = 0;
  Etot_avg[0] = Etot[0];

  // std::ofstream output;
  // int seed = config->time_seed;
  // output.open("errorCheck-"+std::to_string(seed)+".bin", std::ios::binary);

  // for (int i = 0; i < lastStep; i++) {
  //   output.write((char*)&time[i], sizeof(double));
  //   output.write((char*)&work[i], sizeof(double));
  //   output.write((char*)&dE[i], sizeof(double));
  //   output.write((char*)&Etot_avg[i], sizeof(double));
  // }
  // output.close();

  delete [] mass;
  delete [] Etot;
}

// input: cConstants - access time_seed to derive file name
//        generation - for storing into the file, the generation that the mutation is occuring
//        annealing - for storing into the file, the anneal value being used
//        numGenes - the number of genes that are to be mutated
//        recordLog[OPTIM_VARS] - a "mutate mask" that is an array corresponding to the optimization array that contains the random mutate values being applied
// output: file mutateFile-[time_seed].csv is appended a new row containing mutate information
void recordMutateFile(const cudaConstants * cConstants, double generation, double annealing, int numGenes, double recordLog[OPTIM_VARS]) {
  std::ofstream mutateFile;
  int seed = cConstants->time_seed;
  mutateFile.open("mutateFile-" + std::to_string(seed) + ".csv", std::ios_base::app);

  // Record generation, annealing, and number of genes that should be impacted
  mutateFile << generation << "," << annealing << "," << numGenes << ",";

  // Record gamma mutation values
  for (int i = GAMMA_OFFSET; i < (GAMMA_OFFSET + GAMMA_ARRAY_SIZE); i++) {
      mutateFile << recordLog[i] << ",";
  }
  // Record tau mutation values
  for (int i = TAU_OFFSET; i < (TAU_OFFSET + TAU_ARRAY_SIZE); i++) {
      mutateFile << recordLog[i] << ",";
  }
  // Record coast mutation values
  for (int i = COAST_OFFSET; i < (COAST_OFFSET + COAST_ARRAY_SIZE); i++) {
      mutateFile << recordLog[i] << ",";
  }
  // Record alpha, beta, zeta, tripTime
  mutateFile << recordLog[ALPHA_OFFSET] << "," << recordLog[BETA_OFFSET] << "," << recordLog[ZETA_OFFSET] << "," << recordLog[TRIPTIME_OFFSET] << ",";
  mutateFile << "\n";
  
  mutateFile.close();
}

// Output final results of genetic algorithm
// input: x[] - array that holds parameter values of OPTIM_VARS length
//        lastStep - Used to store number of total steps as output
//        threadRank - Rank of the individual element being recorded, currently (July )
//        yOut - Used to store element informatio as output
//        thrust - Passed into rk4sys
//        cConstants - Access constants info such as target element, earth element, derive spaceCraft element, also other values such as rk_tol
// output: yOut contains final eleement information of the spacecraft
//         lastStep contains value last value for number of steps taken
void trajectoryPrint( double x[], int generation, thruster<double> thrust, const cudaConstants* cConstants) {
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
  deltaT = (timeFinal-timeInitial) / cConstants->GuessMaxPossibleSteps; // initial guess for time step, small is preferable

  // setup of thrust angle calculations based off of optimized coefficients
  coefficients<double> coeff;
  initCoefficient(x,coeff, cConstants);
  // Assigning coast threshold (now done in coefficients because is a constant)

  // Assigning wetMass
  double wetMass = cConstants->wet_mass;
  // setting Runge-Kutta tolerance
  double absTol = cConstants->rk_tol;
  // set optmization minimum
  // double Fmin = cConstants->f_min;

  // Initialize memory for the solution vector of the dependant solution
  elements<double>* yp;
  yp = new elements<double>[numSteps];
  
  double *times, *gamma, *tau, *accel_output, *fuelSpent, *work, *dE, *Etot_avg;
  times = new double[numSteps]; // Initialize memory for time array
  gamma = new double[numSteps]; // Initialize memory for gamma array
  tau = new double[numSteps]; // Initialize memory for tau array
  accel_output = new double[numSteps]; // Initialize memory for acceleration array
  fuelSpent = new double[numSteps];  // Initialize memory for fuelSpent array
  work = new double[numSteps];
  dE = new double[numSteps];
  Etot_avg = new double[numSteps];

  double accel; // Initialize memory for  acceleration

  // used to track the cost function throughout a run via output and outputs to a binary
  int lastStepInt;

  rk4sys(timeInitial, x[TRIPTIME_OFFSET] , times, spaceCraft, deltaT, yp, absTol, coeff, accel, gamma, tau, lastStepInt, accel_output, fuelSpent, wetMass, thrust, cConstants);

  double lastStep = lastStepInt;

  // gets the final y values of the spacecrafts for the cost function.
  elements<double> yOut = yp[lastStepInt];

  errorCheck(times, yp, gamma, tau, lastStepInt, accel_output, fuelSpent, wetMass, work, dE, Etot_avg, cConstants);
  progressiveAnalysis(generation, lastStepInt, x, yOut, cConstants);

  std::ofstream output;
  int seed = cConstants->time_seed;
  output.open("orbitalMotion-"+std::to_string(seed)+".bin", std::ios::binary);
  // output.open("orbitalMotion-"+std::to_string(static_cast<int>(seed))+"-"+std::to_string(threadRank)+".bin", std::ios::binary);
  for(int i = 0; i <= lastStepInt; i++) {
    //output << yp[i];
    output.write((char*)&yp[i], sizeof (elements<double>));
    output.write((char*)&times[i], sizeof (double));
    output.write((char*)&gamma[i], sizeof (double));
    output.write((char*)&tau[i], sizeof (double));
    output.write((char*)&accel_output[i], sizeof (double));
    output.write((char*)&fuelSpent[i], sizeof (double));
    output.write((char*)&work[i], sizeof(double));
    output.write((char*)&dE[i], sizeof(double));
    output.write((char*)&Etot_avg[i], sizeof(double));
  }
  output.close();


  double gsize = GAMMA_ARRAY_SIZE, tsize = TAU_ARRAY_SIZE, csize = COAST_ARRAY_SIZE;
  output.open("finalOptimization-"+std::to_string(seed)+".bin", std::ios::binary);
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
  output.write((char*)&earth.r, sizeof(double));
  output.write((char*)&earth.theta, sizeof(double));
  output.write((char*)&earth.z, sizeof(double));
  output.write((char*)&cConstants->fuel_mass, sizeof(double));
  output.write((char*)&cConstants->coast_threshold, sizeof(double));
  output.write((char*)&gsize, sizeof(double));
  output.write((char*)&tsize, sizeof(double));
  output.write((char*)&csize, sizeof(double));

  for (int j = 0; j < OPTIM_VARS; j++) {
    output.write((char*)&x[j], sizeof (double));
  }
  
  output.write((char*)&lastStep, sizeof (double));

  output.close();

  // std::cout << "\ntrajectoryPrint() returned a best posDiff of " << sqrt(pow(cConstants->r_fin_ast - yOut.r, 2) + pow(cConstants->theta_fin_ast - fmod(yOut.theta, 2 * M_PI), 2) + pow(cConstants->z_fin_ast - yOut.z, 2)) << std::endl;
  
  // cleaning up dynamic yp, time, gamma, and tau.
  delete [] yp;
  delete [] times;
  delete [] gamma;
  delete [] tau;
  delete [] accel_output;
  delete [] fuelSpent;
  delete [] work;
  delete [] dE;
  delete [] Etot_avg;
}

// Record progress of individual
// input: output - the output file stream being used
//        rank - the positional performance of the individual
//        ind - the individual object being recorded
//        config - cudaConstants object for accessing thruster_type information
// output: output file is appended information on rank, individual values/parameter information
void progressiveAnalysis(int generation, int numStep, double *start, elements<double> & yp, const cudaConstants *config) {
  int seed = config->time_seed;
  std::ofstream output;
  output.open("progressiveAnalysis.csv", std::ios_base::app);
  output << "\ntime_seed,numStep,posDiff,velDiff,Triptime,alpha,beta,zeta,";

  for (int i = 0; i < GAMMA_ARRAY_SIZE; i++) {
    output << "gamma" << i << ","; 
  }
  for (int i = 0; i < TAU_ARRAY_SIZE; i++) {
    output << "tau" << i << ","; 
  }
  for (int i = 0; i < COAST_ARRAY_SIZE; i++) {
    output << "coast" << i << ","; 
  }
  output << "\n";  

  output << seed << ',' << numStep << ','; 
  output << sqrt(pow(config->r_fin_ast - yp.r, 2) + pow(config->r_fin_ast * config->theta_fin_ast - yp.r * fmod(yp.theta, 2 * M_PI), 2) + pow(config->z_fin_ast - yp.z, 2)) << ',';
  output << sqrt(pow(config->vr_fin_ast - yp.vr, 2) + pow(config->vtheta_fin_ast - yp.vtheta, 2) + pow(config->vz_fin_ast - yp.vz, 2)) << ',';
  output << start[TRIPTIME_OFFSET] << ',' << start[ALPHA_OFFSET] << ',' << start[BETA_OFFSET] << ',' << start[ZETA_OFFSET] << ',';

  for (int i = GAMMA_OFFSET; i < GAMMA_ARRAY_SIZE + GAMMA_OFFSET; i++) {
    output << start[i] << ","; 
  }
  for (int i = TAU_OFFSET; i < TAU_ARRAY_SIZE + TAU_OFFSET; i++) {
    output << start[i] << ","; 
  }
  for (int i = COAST_OFFSET; i < COAST_ARRAY_SIZE + COAST_OFFSET; i++) {
    output << start[i] << ","; 
  }

  output << std::endl;
  output.close();
}

// Initialize the .csv file with header row
// input: cConstants - to access time_seed for deriving file name conventions and also thruster type
// output: file genPerfmanceT-[time_seed].csv, is appended with initial header row info
void initializeRecord(const cudaConstants * cConstants) {
  std::ofstream excelFile;
  int seed = cConstants->time_seed;
  // setting the numeric id tag as the randomization seed (when doing runs of various properties, suggested to add other values to differentiate)
  std::string fileId = std::to_string(seed);
  excelFile.open("genPerformanceT-" + fileId + ".csv", std::ios_base::app);

  excelFile << "Gen,bestPosDiff,bestVelDiff,alpha,beta,zeta,triptime,";
  
  for (int i = 0; i < GAMMA_ARRAY_SIZE; i++) {
    excelFile << "gamma" << i << ","; 
  }
  for (int i = 0; i < TAU_ARRAY_SIZE; i++) {
    excelFile << "tau" << i << ","; 
  }
  for (int i = 0; i < COAST_ARRAY_SIZE; i++) {
    excelFile << "coast" << i << ","; 
  }
  
  
  /*
  for (int i = 0; i < GAMMA_ARRAY_SIZE; i++) {
    excelFile << "gamma" << i << " deviation,";
  }
  for (int i = 0; i < TAU_ARRAY_SIZE; i++) {
    excelFile << "tau" << i << " deviation,";
  }
  excelFile << "alpha deviation, beta deviation, zeta deviation, triptime deviation,";
  
  for (int i = 0; i < COAST_ARRAY_SIZE; i++) {
    excelFile << "coast" << i << " deviation,";
  }*/

  excelFile << ",\n";
  excelFile.close();
    // call setMutateFile to set it up
    // setMutateFile(cConstants);
}

// Take in the current state of the generation and appends to files
// input: cConstants - access time_seed to derive file name
//        pool - passes pool[0] to writeIndividualToFiles & writeThrustToFiles
//        generation - passed to writeIndividualToFiles & writeThrustToFiles
//        new_anneal - passed into writeIndividualToFiles
//        poolSize - used to access worst individual in pool by using index value poolSize-1
//        thrust - used in conditional statement of thrust type
// output: files BestInGenerations/WorstInGenerations have appended information using writeIndividualToFiles method
//         files BestThrustGens/WorstThurstGens have appended information using writeThrustFiles method
void recordGenerationPerformance(const cudaConstants * cConstants, Individual * pool, double generation, double new_anneal, int poolSize, thruster<double>& thrust) {
  std::ofstream excelFile;
  int seed = cConstants->time_seed;
  std::string fileId = std::to_string(seed);
  excelFile.open("genPerformanceT-" + fileId + ".csv", std::ios_base::app);

  excelFile << generation << "," << pool[0].posDiff << ",";
  excelFile << pool[0].velDiff << ",";

  excelFile << pool[0].startParams.alpha << ",";
  excelFile << pool[0].startParams.beta << ",";
  excelFile << pool[0].startParams.zeta << ",";
  excelFile << pool[0].startParams.tripTime << ",";

  for (int i = GAMMA_OFFSET; i < GAMMA_ARRAY_SIZE + GAMMA_OFFSET; i++) {
    excelFile << pool[0].startParams.coeff.gamma[i-GAMMA_OFFSET] << ","; 
  }
  for (int i = TAU_OFFSET; i < TAU_ARRAY_SIZE + TAU_OFFSET; i++) {
    excelFile << pool[0].startParams.coeff.tau[i-TAU_OFFSET] << ","; 
  }
  for (int i = COAST_OFFSET; i < COAST_ARRAY_SIZE + COAST_OFFSET; i++) {
    excelFile << pool[0].startParams.coeff.coast[i-COAST_OFFSET] << ","; 
  }

/*
  double * deviation_array = new double[OPTIM_VARS];
  generationDiversity(pool, poolSize, deviation_array);
  
  for (int gene = 0; gene < OPTIM_VARS; gene++) {
    excelFile << deviation_array[gene] << ",";
  }
  delete [] deviation_array;*/
  excelFile << "\n"; // End of row
  excelFile.close();
}


// Takes in a pool and records the parameter info on all individuals
// input: cConstants - to access time_seed in deriving file name
//        pool - holds all the individuals to be stored
//        poolSize - to use in iterating through the pool
//        generation - used in deriving file name
// output: file generation#[generation]-[time_seed].csv is created with each row holding parameter values of individuals
void recordAllIndividuals(const cudaConstants * cConstants, Individual * pool, int poolSize, int generation) {
  std::ofstream entirePool;
  int seed = cConstants->time_seed;
  entirePool.open("generation#" + std::to_string(generation) + "-" + std::to_string(seed) + ".csv");
  // Setup the header row
  entirePool << "position,alpha,beta,zeta,tripTime,";
  for (int i = 0; i < GAMMA_ARRAY_SIZE; i++) {
    entirePool << "gamma" << i << ",";
  }
  for (int i = 0; i < TAU_ARRAY_SIZE; i++) {
    entirePool << "tau" << i << ",";
  }
  for (int i = 0; i < COAST_ARRAY_SIZE; i++) {
    entirePool << "coast" << i << ",";
  }
  entirePool << '\n';

  // Record all individuals in the pool
  for (int i = 0; i < poolSize; i++) {
    entirePool << i << ",";
    entirePool << pool[i].startParams.alpha << ",";
    entirePool << pool[i].startParams.beta << ",";
    entirePool << pool[i].startParams.zeta << ",";
    entirePool << pool[i].startParams.tripTime << ",";

    for (int j = 0; j < GAMMA_ARRAY_SIZE; j++) {
      entirePool << pool[i].startParams.coeff.gamma[j] << ",";
    }
    for (int j = 0; j < TAU_ARRAY_SIZE; j++) {
      entirePool << pool[i].startParams.coeff.tau[j] << ",";
    }
    for (int j = 0; j < COAST_ARRAY_SIZE; j++) {
      entirePool << pool[i].startParams.coeff.coast[j] << ",";
    }
    entirePool << "\n";
  }
  entirePool.close();
}


// Method for doing recording information at the end of the optimization process
// input: cConstants - access record_mode, if record_mode == true then call progressiveRecord method, also passed into writeTrajectoryToFile method as well as progressiveRecord
//        pool - To access the best individual (pool[0])
//        generation - to record the generation value 
//        thrust - passed into progressiveRecord and writeTrajectoryToFile
// output: writeTrajectoryToFile is called, if in record_mode then progressiveRecord is called as well
void finalRecord(const cudaConstants* cConstants, Individual * pool, int generation, thruster<double>& thrust) {
  // To store parameter values and pass onto writeTrajectoryToFile
  double *start = new double[OPTIM_VARS];

  // Output the final best individual
  for (int j = 0; j < pool[0].startParams.coeff.gammaSize; j++) {
      start[GAMMA_OFFSET + j] = pool[0].startParams.coeff.gamma[j];
  }
  for (int j = 0; j < pool[0].startParams.coeff.tauSize; j++) {
      start[TAU_OFFSET + j] = pool[0].startParams.coeff.tau[j];
  }
  for (int j = 0; j < pool[0].startParams.coeff.coastSize; j++) {
      start[COAST_OFFSET + j] = pool[0].startParams.coeff.coast[j];
  }

  start[TRIPTIME_OFFSET] = pool[0].startParams.tripTime;
  start[ALPHA_OFFSET] = pool[0].startParams.alpha;
  start[BETA_OFFSET] = pool[0].startParams.beta;
  start[ZETA_OFFSET] = pool[0].startParams.zeta;

  std::cout << "Comparison\n";
  std::cout << "CUDA posDiff: " << pool[0].posDiff << std::endl;
  std::cout << "CUDA velDiff: " << pool[0].velDiff << std::endl;

  // Could instead use a ratio between position and velocity differnce as done in comparison of Individuals
  trajectoryPrint(start, generation, thrust, cConstants);

  // std::cout << "\nfinalRecord() returned a best posDiff of " << pool[0].posDiff << std::endl;

  delete [] start;
}

// method that stores information of launchCon of timeRes*24 resolution
// input: cConstants - access time range and resolution info on launchCon
//        launchCon - access elements of earth 
// output: EarthCheckValues.csv is created and holds rows of element info on earth with timeStamp on each row
void recordEarthData(const cudaConstants * cConstants) {
  double timeStamp = cConstants->triptime_min*SECONDS_IN_YEAR;

  int seed = cConstants->time_seed;
  std::ofstream earthValues;
  earthValues.open("EarthCheckValues-"+ std::to_string(seed) +".csv");
  // Set header row for the table to record values, with timeStamp
  earthValues << "TimeStamp, Radius, Theta, Z, vRadius, vTheta, vZ\n";
  while (timeStamp < cConstants->triptime_max*SECONDS_IN_YEAR) {
      earthValues << timeStamp << "," << launchCon->getCondition(timeStamp);
      timeStamp += cConstants->timeRes*24; // Increment to next day as timeRes is assumed to be set to every hour in this output
  }
  // Done recording earth calculations, close file
  earthValues.close();
}