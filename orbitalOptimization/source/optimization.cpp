//Didymos-Optimization_Project:
//Last Editor: Ben and Mateo
//Tasks Completed: 
  //Created iterativeOptimization() and optimizeStartConditions().
  //Changed all starts and steps to be defined constants instead of magic numbers.
  //Changed the initial guess of a parameter to be a random number within a resonable range of values.


#include "optimization.h" 
#include "nelder_mead.h" // used for nelmin()
#include "constants.h" //used for wetMass
#include "orbitalMotion.h" //used for trajectory() and trajectoryPrint()
#include <iostream> // cout
#include <iomanip> //used for setw(), sets spaces between values
#include <time.h> //for seeding the random number generator
#include <random>
#include <chrono>


int main ()
{
  
  //iterativeOptimize(); // manually set initial conditions
  optimizeStartConditions(); // random values within a given range for initial conditions

  return 0;
}

void optimizeStartConditions(){
  //Pre-allocating memory for starting parameters and steps
  //start - starting parameters for optimization, order and contents defined in constants.h
  //step - starting step sizes for starting parameters
  double *start = new double[OPTIM_VARS];
  double *step = new double[OPTIM_VARS];

  double bestCost = 1.0E9;
  //double *bestStart = new double[OPTIM_VARS];

  std::mt19937 mt_rand(time(0)); //seed the random number generator

  std::ofstream output;
  output.open ("optimized-start-conditions.txt");

  int executions = 80;
  for(int i = 0; i < executions; i++)
  {
    // Initial guesses for variables based off of previous runs which have small cost values
    start[GAMMA_OFFSET] = mt_rand() % 201/10.0 - 10.0; // -10 - 10
    start[GAMMA_OFFSET+1] = mt_rand() % 201/10.0 - 10.0;
    start[GAMMA_OFFSET+2] = mt_rand() % 201/10.0 - 10.0;
    start[GAMMA_OFFSET+3] = mt_rand() % 201/10.0 - 10.0;
    start[GAMMA_OFFSET+4] = mt_rand() % 201/10.0 - 10.0;
    start[GAMMA_OFFSET+5] = mt_rand() % 201/10.0 - 10.0;
    start[GAMMA_OFFSET+6] = mt_rand() % 201/10.0 - 10.0;

    start[TAU_OFFSET] = mt_rand() % 201/10.0 - 10.0; // -10.0 - 10.0
    start[TAU_OFFSET+1] = mt_rand() % 201/10.0 - 10.0;
    start[TAU_OFFSET+2] = mt_rand() % 201/10.0 - 10.0;

    start[ALPHA_OFFSET] = (mt_rand() % 629) / 100.0 - 3.14; // -pi - pi
    start[BETA_OFFSET] = (mt_rand() % 629) / 100.0 - 3.14;
    start[ZETA_OFFSET] = (mt_rand() % 315) / 100.0 - 1.57;

    start[TRIPTIME_OFFSET] = 365*24*3600*(std::rand() % 10001 / 10000.0 + 1.0); // 1.0 - 2.0 years converted to seconds

    start[COAST_OFFSET] = mt_rand() % 201/10.0 - 10.0; // -10.0 - 10.0
    start[COAST_OFFSET+1] = mt_rand() % 201/10.0 - 10.0;
    start[COAST_OFFSET+2] = mt_rand() % 201/10.0 - 10.0;
    start[COAST_OFFSET+3] = mt_rand() % 201/10.0 - 10.0;
    start[COAST_OFFSET+4] = mt_rand() % 201/10.0 - 10.0;

    // Initial change in variable size based on the variable start value
    // Delimits the search space
    step[GAMMA_OFFSET] = 1.0E01/2;
    step[GAMMA_OFFSET+1] = 1.0E01/2;
    step[GAMMA_OFFSET+2] = 1.0E01/2;
    step[GAMMA_OFFSET+3] = 1.0E01/2;
    step[GAMMA_OFFSET+4] = 1.0E01/2;
    step[GAMMA_OFFSET+5] = 1.0E01/2;
    step[GAMMA_OFFSET+6] = 1.0E01/2;

    step[TAU_OFFSET] = 1.0E0;
    step[TAU_OFFSET+1] = 1.0E0;
    step[TAU_OFFSET+2] = 1.0E0;

    step[ALPHA_OFFSET] = 1.0E00;
    step[BETA_OFFSET] = 1.0E00;
    step[ZETA_OFFSET] = 1.0E00;

    step[TRIPTIME_OFFSET] = 1.0E07;

    step[COAST_OFFSET] = 1.0E01;
    step[COAST_OFFSET+1] = 1.0E01;
    step[COAST_OFFSET+2] = 1.0E01;
    step[COAST_OFFSET+3] = 1.0E01;
    step[COAST_OFFSET+4] = 1.0E01;


    optimizing(start, step);
    std::cout<<"No problems in optimizing\n";
    // writes the solution based on optimized variables to a binary file
    int numSteps = 0;
    double cost; // to store the cost caluclated by trajectoryPrint()

    writeTrajectoryToFile(start, cost);

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
      // code not outputing the right start values
      //bestStart = start;
    }

    std::cout << "run " << i + 1 << " complete" << std::endl;
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

  std::mt19937 mt_rand(time(0)); //seed the random number generator

  // Initial guesses for variables based off of previous runs which have small cost values
      // Initial guesses for variables based off of previous runs which have small cost values
  start[GAMMA_OFFSET] = mt_rand() % 201/10.0 - 10.0; // -10 - 10
  start[GAMMA_OFFSET+1] = mt_rand() % 201/10.0 - 10.0;
  start[GAMMA_OFFSET+2] = mt_rand() % 201/10.0 - 10.0;
  start[GAMMA_OFFSET+3] = mt_rand() % 201/10.0 - 10.0;
  start[GAMMA_OFFSET+4] = mt_rand() % 201/10.0 - 10.0;
  start[GAMMA_OFFSET+5] = mt_rand() % 201/10.0 - 10.0;
  start[GAMMA_OFFSET+6] = mt_rand() % 201/10.0 - 10.0;


  start[TAU_OFFSET] = mt_rand() % 201/10.0 - 10.0; // -10.0 - 10.0
  start[TAU_OFFSET+1] = mt_rand() % 201/10.0 - 10.0;
  start[TAU_OFFSET+2] = mt_rand() % 201/10.0 - 10.0;

  start[ALPHA_OFFSET] = (mt_rand() % 629) / 100.0 - 3.14; // -pi - pi
  start[BETA_OFFSET] = (mt_rand() % 629) / 100.0 - 3.14;
  start[ZETA_OFFSET] = (mt_rand() % 315) / 100.0 - 1.57;

  start[TRIPTIME_OFFSET] = 365*24*3600*(std::rand() % 10001 / 10000.0 + 1); // 1.5 - 2.5 years converted to seconds

  start[COAST_OFFSET] = mt_rand() % 201/10.0 - 10.0; // -10.0 - 10.0
  start[COAST_OFFSET+1] = mt_rand() % 201/10.0 - 10.0;
  start[COAST_OFFSET+2] = mt_rand() % 201/10.0 - 10.0;
  start[COAST_OFFSET+3] = mt_rand() % 201/10.0 - 10.0;
  start[COAST_OFFSET+4] = mt_rand() % 201/10.0 - 10.0;

  // Initial change in variable size based on the variable start value
  // Delimits the search space
  step[GAMMA_OFFSET] = 1.0E01/2;
  step[GAMMA_OFFSET+1] = 1.0E01/2;
  step[GAMMA_OFFSET+2] = 1.0E01/2;
  step[GAMMA_OFFSET+3] = 1.0E01/2;
  step[GAMMA_OFFSET+4] = 1.0E01/2;
  step[GAMMA_OFFSET+5] = 1.0E01/2;
  step[GAMMA_OFFSET+6] = 1.0E01/2;

  step[TAU_OFFSET] = 1.0E0;
  step[TAU_OFFSET+1] = 1.0E0;
  step[TAU_OFFSET+2] = 1.0E0;

  step[ALPHA_OFFSET] = 1.0E00;
  step[BETA_OFFSET] = 1.0E00;
  step[ZETA_OFFSET] = 1.0E00;

  step[TRIPTIME_OFFSET] = 1.0E07;

  step[COAST_OFFSET] = 1.0E01;
  step[COAST_OFFSET+1] = 1.0E01;
  step[COAST_OFFSET+2] = 1.0E01;
  step[COAST_OFFSET+3] = 1.0E01;
  step[COAST_OFFSET+4] = 1.0E01;

  // For loop to reutilize the final value of the c vector as the guess for the next optimization 
  int executions = 10;
  for(int i = 0; i < executions; i++)
  {
    optimizing(start, step);
  }

  // writes the solution based on optimized variables to a binary file
  double cost = 0;
  writeTrajectoryToFile(start,cost);

  delete [] start;
  delete [] step;
}

void optimizing (double *&start, double *step)
//  Purpose: Optimize the following:
//* Coefficients for gamma and tau fourier series,
//* alpha and beta angles (used in initial velocity of the spacecraft),
//* trip times,
//* coast fourier
//* coast threshold
{
  // initializing variables for nelmin algorithm. See nelder_mead.cpp for input/output information
  int i; // iteration number
  int icount; // number of function evaluations
  int ifault; // Output, int *IFAULT, error indicator. 0, no errors detected. 1, REQMIN, N, or KONVGE has an illegal value. 
  //2, iteration terminated because KCOUNT was exceeded without convergence.
  int kcount; // maximum number of iterations
  int konvge; // how often cost value is checked for convergence
  int numres; // number of restarts
  double reqmin; // the convergence minimum
  double *xmin; // final output
  double ynewlo; // cost value

  //Allocting xmin: space for nelder_mead algorithm to fill in final optimized parameters
  xmin = new double[OPTIM_VARS];

  // Terminating limit for the variance of function values
  // nelmin algorithm aims for the square root of this number
  reqmin = 1.0E-40; // Expecting a solution with cost within 10E-20 error
  
  // how often the equation checks for a convergence
  konvge = 20+std::rand()%2;
  // maximum number of iterations for convergence
  kcount = 10000+std::rand()%100;

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


void writeTrajectoryToFile(double *start, double & cost)
{
  int numSteps = 0;

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
}