#include "optimization.h" 
#include "nelder_mead.h" // used for nelmin()
#include "constants.h" //used for wetMass
#include "orbitalMotion.h" //used for trajectory() and trajectoryPrint()
#include <iostream> // cout
#include <iomanip> //used for setw(), sets spaces between values
#include <time.h> //for seeding the random number generator


int main ()
{
  timestamp ();

  //optimizing();
  //iterativeOptimize();
  optimizeStartConditions();

  timestamp ();
  
  return 0;
}

void optimizeStartConditions(){
  // allocating memory according to number of variables
  double *start = new double[OPTIM_VARS];
  double *step = new double[OPTIM_VARS];

  double bestCost = 1.0E9;
  //double *bestStart = new double[OPTIM_VARS];

  std::srand(std::time(NULL)); //seed the random number generator

  std::ofstream output;
  output.open ("optimized-start-conditions.txt");

  for(int i = 0; i < 10000; i++){
    // random initial guesses for variables within a reasonable range
    start[GAMMA_OFFSET] = std::rand() % 201 - 100; // -100 - 100
    start[GAMMA_OFFSET+1] = std::rand() % 201 - 100;
    start[GAMMA_OFFSET+2] = std::rand() % 201 - 100;
    start[GAMMA_OFFSET+3] = std::rand() % 201 - 100;
    start[GAMMA_OFFSET+4] = std::rand() % 201 - 100;
    start[GAMMA_OFFSET+5] = std::rand() % 201 - 100;
    start[GAMMA_OFFSET+6] = std::rand() % 201 - 100;
    start[GAMMA_OFFSET+7] = std::rand() % 201 - 100;
    start[GAMMA_OFFSET+8] = std::rand() % 201 - 100;

    start[TAU_OFFSET] = (std::rand() % 201) / 10 - 10; // -10.0 - 10.0
    start[TAU_OFFSET+1] = (std::rand() % 201) / 10 - 10;
    start[TAU_OFFSET+2] = (std::rand() % 201) / 10 - 10;

    start[ALPHA_OFFSET] = (std::rand() % 201) / 10 - 10;
    start[BETA_OFFSET] = (std::rand() % 201) / 10 - 10;

    start[TRIPTIME_OFFSET] = 365*24*3600*(std::rand() % 50001 / 10000 + 1); // 1.0000 - 6.0000 years converted to seconds

    start[COAST_OFFSET] = (std::rand() % 101) / 10; //0.0 - 10.0
    start[COAST_OFFSET+1] = (std::rand() % 101) / 10;
    start[COAST_OFFSET+2] = (std::rand() % 101) / 10;
    start[COAST_OFFSET+3] = (std::rand() % 101) / 10;
    start[COAST_OFFSET+4] = (std::rand() % 101) / 10;

    start[THRESHOLD_OFFSET] = (std::rand() % 101) / 10;
    //start[WETMASS_OFFSET] = DRY_MASS+200; // 3950 kg

    // initial change in variable size
    // based on the variable start value
    step[0] = 1.0E01;
    step[1] = 1.0E01;
    step[2] = 1.0E01;
    step[3] = 1.0E01;
    step[4] = 1.0E01;
    step[5] = 1.0E01;
    step[6] = 1.0E01;
    step[7] = 1.0E01;
    step[8] = 1.0E01;
    step[9] = 1.0E01;
    step[10] = 1.0E01;
    step[11] = 1.0E01;
    step[12] = 1.0E00;
    step[13] = 1.0E00;
    step[14] = 1.0E07;
    step[15] = 1.0E00;
    step[16] = 1.0E00;
    step[17] = 1.0E00;
    step[18] = 1.0E00;
    step[19] = 1.0E00;
    step[20] = 1.0E-02;
    //step[21] = 1.0E01;

    optimizing(start, step);

    // writes the solution based on optimized variables to a binary file
    int numSteps = 0;
    double cost; // to store the cost caluclated by trajectoryPrint()

    trajectoryPrint(start, numSteps, cost);

    output << "start values:" << std::endl;
    for(int i = 0; i < OPTIM_VARS / 2 + 1; i++)
    {
      output << i + 1 << ": " << start[i] << ", ";
    }
    output << std::endl;
    for(int i = OPTIM_VARS / 2 + 1; i < OPTIM_VARS; i++)
    {
      output << i + 1<< ": " << start[i] << ", ";
    }
    output << std::endl << "cost value: " << cost << std::endl;
    output << "---------------------------------------------------------------------------------" << std::endl;

    if(cost < bestCost){
      bestCost = cost;
      // not outputing the right start values
      //bestStart = start;
    }
  }
  output << "---------------------------------------------------------------------------------" << std::endl;
  output << "---------------------------------------------------------------------------------" << std::endl;
  output << "BEST RESULTS:" << std::endl;
  /*
  output << "start values:" << std::endl;
  for(int i = 0; i < OPTIM_VARS / 2 + 1; i++)
  {
    output << i + 1 << ": " << bestStart[i] << ", ";
  }
  output << std::endl;
  for(int i = OPTIM_VARS / 2 + 1; i < OPTIM_VARS; i++)
  {
    output << i + 1<< ": " << bestStart[i] << ", ";
  }
  */
  output << std::endl << "cost value: " << bestCost << std::endl;
  output.close();

  delete [] start;
  //delete [] bestStart;
  delete [] step;
}

void iterativeOptimize(){
  // allocating memory according to number of variables
  double *start = new double[OPTIM_VARS];
  double *step = new double[OPTIM_VARS];

  // x[0]-x[8]: gamma coefficients used to calculate fourier series
  // x[9]-x[11]: tau coefficients used to calculate fourier series
  // x[12]: alpha - launch angle (declination) position 
  // x[13]: beta - launch angle (declination) velocity 
  // x[14]: trip time - total time from launch to impact, sets the initial earth position
  // x[15-19]: coast coefficients used to calculate fourier series
  // x[20]: coast threshold - value set to determine when coasting occurs
  // x[21]: wet mass - mass of spacecraft including fuel
  // Initial guesses for variables based off of previous runs which have small cost values
  start[GAMMA_OFFSET] = 10;
  start[GAMMA_OFFSET+1] = 10;
  start[GAMMA_OFFSET+2] = 10;
  start[GAMMA_OFFSET+3] = 10;
  start[GAMMA_OFFSET+4] = 10;
  start[GAMMA_OFFSET+5] = 10;
  start[GAMMA_OFFSET+6] = 10;
  start[GAMMA_OFFSET+7] = 10;
  start[GAMMA_OFFSET+8] = 10;
  start[TAU_OFFSET] = 10;
  start[TAU_OFFSET+1] = 10;
  start[TAU_OFFSET+2] = 10;
  start[TAU_OFFSET+3] = 10;
  start[TAU_OFFSET+4] = 10;
  start[ALPHA_OFFSET] = 0.5;
  start[BETA_OFFSET] = 0.5;
  start[TRIPTIME_OFFSET] = 365*24*3600*1.5; // 2 YEARS
  start[COAST_OFFSET] = 0.5;
  start[COAST_OFFSET+1] = 0.5;
  start[COAST_OFFSET+2] = 0.5;
  start[COAST_OFFSET+3] = 0.5;
  start[COAST_OFFSET+4] = 0.5;
  start[THRESHOLD_OFFSET] = 0.05;
  //start[WETMASS_OFFSET] = DRY_MASS+200; // 3950 kg

  // Initial change in variable size based on the variable start value
  // Delimits the search space
  step[GAMMA_OFFSET] = 1.0E02;
  step[GAMMA_OFFSET+1] = 1.0E02;
  step[GAMMA_OFFSET+2] = 1.0E02;
  step[GAMMA_OFFSET+3] = 1.0E02;
  step[GAMMA_OFFSET+4] = 1.0E02;
  step[GAMMA_OFFSET+5] = 1.0E02;
  step[GAMMA_OFFSET+6] = 1.0E02;
  step[GAMMA_OFFSET+7] = 1.0E02;
  step[GAMMA_OFFSET+8] = 1.0E02;
  step[TAU_OFFSET] = 1.0E02;
  step[TAU_OFFSET+1] = 1.0E02;
  step[TAU_OFFSET+2] = 1.0E02;
  step[TAU_OFFSET+3] = 1.0E02;
  step[TAU_OFFSET+4] = 1.0E02;
  step[ALPHA_OFFSET] = 1.0E00;
  step[BETA_OFFSET] = 1.0E00;
  step[TRIPTIME_OFFSET] = 1.0E07;
  step[COAST_OFFSET] = 1.0E02;
  step[COAST_OFFSET+1] = 1.0E02;
  step[COAST_OFFSET+2] = 1.0E02;
  step[COAST_OFFSET+3] = 1.0E02;
  step[COAST_OFFSET+4] = 1.0E02;
  step[THRESHOLD_OFFSET] = 1.0E-02;
  //step[21] = 1.0E01;

  // For loop to reutilize the final value of the c vector as the guess for the next optimization 
  for(int i = 0; i < 1; i++)
  {
    optimizing(start, step);
  }

 // writes the solution based on optimized variables to a binary file
  int numSteps = 0;
  double cost = 0;

  trajectoryPrint(start, numSteps, cost);

  //writes final optimization values to a seperate file
  std::ofstream output;

  output.open ("final-optimization.bin", std::ios::binary);
  for(int i=0; i < OPTIM_VARS; i++)
  {
    output.write((char*)&start[i], sizeof (double));
  }
  output.write((char*)&numSteps, sizeof (int));
  output.close();

  delete [] start;
  delete [] step;
}

void optimizing (double *&start, double *step)
//  Purpose: Optimize the following:
//* Coefficients for gamma and tau fourier series,
//* alpha and beta angles (used in initial velocity of the spacecraft),
//* trip times,
//* Fuel mass.
{
  // initializing variables for nelmin algorithm. See nelder_mead.cpp for input/output information
  int i; 
  int icount; 
  int ifault;
  int kcount;
  int konvge;
  int numres;
  double reqmin;
  double *xmin;
  double ynewlo;

  //Allocting xmin: space for nelder_mead algorithm to fill in final optimized parameters
  xmin = new double[OPTIM_VARS];

  // Terminating limit for the variance of function values
  // nelmin algorithm aims for the square root of this number
  reqmin = 1.0E-40; // Expecting a solution with cost within 10E-20 error
  
  // how often the equation checks for a convergence
  konvge = 20+std::rand()%2;
  // maximum number of iterations for convergence
  kcount = 30000+std::rand()%100;

    //****************
    // Move into its own function
    std::cout << "\n"<<"Starting conditions:"<<std::endl;
    for ( i = 0; i < OPTIM_VARS; i++ )
    {
      std::cout << std::setw(2) << start[i] << ", ";
    }

    // optimization value for the initial conditions
    ynewlo = trajectory (start);

    std::cout << "\n"<< " F(X) = " << ynewlo << std::endl;
  
  // nelder_mead function (optimization function)
  // see nelder_mead.cpp for input and output information
  nelmin (trajectory, OPTIM_VARS, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);

  //****************
  // Move into its own function
  // displays error type when an error occurs
  std::cout << "\nReturn code IFAULT = " << ifault << "\nEstimate of minimizing value X*:\n\n";
  for (i = 0; i < OPTIM_VARS; i++)
  {
    std::cout << std::setw(2) << xmin[i] << ",";
  }
  std::cout << "\nF(X) = " << ynewlo << "\n";
  std::cout << "\n"<< "  Number of iterations = " << icount << "\n"<< "  Number of restarts =   " << numres << "\n";

  // use the results as the starting point for the next run
  delete [] start;
  start = xmin;

  return;
}
