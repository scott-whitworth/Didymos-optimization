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
#include "earthInfo.h"
#include <iostream> // cout
#include <iomanip> //used for setw(), sets spaces between values
#include <time.h> //for seeding the random number generator
#include <random>
#include <chrono>
#include <string>
#include <fstream>


int main ()
{
  //////////////////////////////////////////////////////////////////////////////////
  //Global variable needs to be initialized
  // This is used to get the earth position for each trip time
  // Define variables to be passed into EarthInfo
  double startTime = 15778800; // 0.5 year (s)
  double endTime = 78894000; // 2.5 years (s)
  double timeRes = 3600; // position of earth is calculated for every hour

  // initializes EarthInfo
  launchCon = new EarthInfo(startTime, endTime, timeRes);
  ////////////////////////////////////////////////////////////////////////////////////

  
  int numTries = 2;
  //optimizeStartConditions(numTries); // random values within a given range for initial conditions

  checkBinaryFile(14);

  //iterativeOptimize(); // manually set initial conditions

  delete launchCon;
  return 0;
}

void optimizeStartConditions(int executions)
{
  //Pre-allocating memory for starting parameters and steps
  //start - starting parameters for optimization, order and contents defined in constants.h
  //step - starting step sizes for starting parameters
  double *start = new double[OPTIM_VARS];
  double *step = new double[OPTIM_VARS];

  std::mt19937 mt_rand(time(0)); //seed the random number generator

  std::ofstream output;
  output.open ("optimized-start-conditions.txt");

  for(int i = 0; i < executions; i++)
  {
    // Initial guesses for variables based off of previous runs which have small cost values
    start[GAMMA_OFFSET] = mt_rand() % 201/100.0 - 1.0; // -10 - 10
    start[GAMMA_OFFSET+1] = mt_rand() % 201/100.0 - 1.0;
    start[GAMMA_OFFSET+2] = mt_rand() % 201/100.0 - 1.0;
    start[GAMMA_OFFSET+3] = mt_rand() % 201/100.0 - 1.0;
    start[GAMMA_OFFSET+4] = mt_rand() % 201/100.0 - 1.0; 
    start[GAMMA_OFFSET+5] = mt_rand() % 201/100.0 - 1.0;
    start[GAMMA_OFFSET+6] = mt_rand() % 201/100.0 - 1.0;

    start[TAU_OFFSET] = mt_rand() % 201/100.0 - 1.0; // -10.0 - 10.0
    start[TAU_OFFSET+1] = mt_rand() % 201/100.0 - 1.0;
    start[TAU_OFFSET+2] = mt_rand() % 201/100.0 - 1.0;

    start[ALPHA_OFFSET] = (mt_rand() % 629) / 100.0 - 3.14; // -pi - pi
    start[BETA_OFFSET] = (mt_rand() % 629) / 100.0 - 3.14;
    start[ZETA_OFFSET] = (mt_rand() % 315) / 100.0 - 1.57;

    start[TRIPTIME_OFFSET] = 365*24*3600*(std::rand() % 10001 / 10000.0 + 1.0); // 1.0 - 2.0 years converted to seconds

    start[COAST_OFFSET] = mt_rand() % 201/100.0 - 1.0; // -10.0 - 10.0
    start[COAST_OFFSET+1] = mt_rand() % 201/100.0 - 1.0;
    start[COAST_OFFSET+2] = mt_rand() % 201/100.0 - 1.0;
    start[COAST_OFFSET+3] = mt_rand() % 201/100.0 - 1.0;
    start[COAST_OFFSET+4] = mt_rand() % 201/100.0 - 1.0;

    // Initial change in variable size based on the variable start value
    // Delimits the search space

    setStep(step);

    optimizing(start, step);
    // writes the solution based on optimized variables to a binary file
    double cost = trajectory(start); // to store the cost caluclated by trajectoryPrint()

    if(trajectory(start)<2*pow(10,-20))
    {
      writeTrajectoryToFile(start, cost,i);
    }
    
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
    std::cout << "run " << i + 1 << " complete" << std::endl;
  }
  output << "---------------------------------------------------------------------------------" << std::endl;
  output << "---------------------------------------------------------------------------------" << std::endl;

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
  //start[GAMMA_OFFSET] = 0.273664,-1.0756,0.6822,-1.27997,-1.45467,-0.0136089,-0.884055,-0.162799,0.00384891,1.04031,-0.726513,1.07761,-1.49482,3.27689e+007,-0.736647,-0.763116,1.37091,0.0260613,-0.418265,; // -10 - 10
  start[GAMMA_OFFSET+1] = -2.7824020148443617816980122370e-01;
  start[GAMMA_OFFSET+2] = 3.9706971844867466892026186542e-01;
  start[GAMMA_OFFSET+3] = -2.6246406920124598638466295597e-01;
  start[GAMMA_OFFSET+4] = -1.6050991385640158704006807966e+00;
  start[GAMMA_OFFSET+5] = 3.1371729669213965774332564251e-01;
  start[GAMMA_OFFSET+6] = 3.7885357138619657479949864864e-01;


  start[TAU_OFFSET] = -3.3390650367959540112394734024e-01; // -5.0 - 5.0
  start[TAU_OFFSET+1] = 5.5499614221243376288583704081e-01;
  start[TAU_OFFSET+2] = 5.1361430331364199552979243890e-01;

  start[ALPHA_OFFSET] = 1.7574222525172265019222095361e+00; // -pi - pi
  start[BETA_OFFSET] = -1.6837629692054056906869163868e+00;
  start[ZETA_OFFSET] = 9.8003867179619408300794702882e-01;

  start[TRIPTIME_OFFSET] = 2.8695261326080441474914550781e+07; // 1.5 - 2.5 years converted to seconds

  start[COAST_OFFSET] = -8.6302487703370944771563699760e-01; // -10.0 - 10.0
  start[COAST_OFFSET+1] = 4.1165731650421222287405953466e-01;
  start[COAST_OFFSET+2] = 2.8647589778377002822651320457e-01;
  start[COAST_OFFSET+3] = -4.2160094983718154892926577304e-01;
  start[COAST_OFFSET+4] = 2.5491324082429184239018127300e-01;



  setStep(step);

  // For loop to reutilize the final value of the c vector as the guess for the next optimization 
  int executions = 1;
  for(int i = 0; i < executions; i++)
  {
    optimizing(start, step);
  }

  // writes the solution based on optimized variables to a binary file
  double cost = 0;
  writeTrajectoryToFile(start,cost,0);

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

// Todo write how it works
void writeTrajectoryToFile(double *start, double & cost, int i)
{
  double numStep = 0;
  elements<double> yp;
  trajectoryPrint(start, numStep, cost,i,yp);
  //writes final optimization values to a seperate file
  std::ofstream output;

  output.open ("final-optimization"+std::to_string(i)+".bin", std::ios::binary);
  for(int j=0; j < OPTIM_VARS; j++)
  {
    output.write((char*)&start[j], sizeof (double));
  }
  output.write((char*)&numStep, sizeof (double));
  output.close();
}

// Initial change in variable size based on the variable start value
// Delimits the search space
void setStep(double step[])
{
  step[GAMMA_OFFSET] = 1.0E00;
  step[GAMMA_OFFSET+1] = 1.0E00;
  step[GAMMA_OFFSET+2] = 1.0E00;
  step[GAMMA_OFFSET+3] = 1.0E00;
  step[GAMMA_OFFSET+4] = 1.0E00;
  step[GAMMA_OFFSET+5] = 1.0E00;
  step[GAMMA_OFFSET+6] = 1.0E00;

  step[TAU_OFFSET] = 5.0E-01;
  step[TAU_OFFSET+1] = 5.0E-01;
  step[TAU_OFFSET+2] = 5.0E-01;

  step[ALPHA_OFFSET] = 1.0E00;
  step[BETA_OFFSET] = 1.0E00;
  step[ZETA_OFFSET] = 1.0E00;

  step[TRIPTIME_OFFSET] = 1.0E07;

  step[COAST_OFFSET] = 5.0E-01;
  step[COAST_OFFSET+1] = 5.0E-01;
  step[COAST_OFFSET+2] = 5.0E-01;
  step[COAST_OFFSET+3] = 5.0E-01;
  step[COAST_OFFSET+4] = 5.0E-01;
}

// If you want to check the results saved in a binary file
// Size is just the number of solutions of the file
// The code is also expecting 19 parameters in the parameter vector
void checkBinaryFile(int size)
{
    
  std::ifstream starts;
  starts.open("optimizedVector.bin", std::ifstream::in|std::ios::binary);

  double startDoubles;

  // sort the data into 2 dimensions
  // one row is one set of starting parameters
  // each column is a specific variable:
  // 0-6 gamma
  // 7-9 tau
  // 10-12 launch angles
  // 13 trip time
  // 14-19 coast
  int parameterSize = 19;

  double arrayCPU[size][parameterSize];
  double singleArray[parameterSize];
  for(int i = 0; i < parameterSize; i++)
  {  // rows
    for(int j = 0; j < size; j++)
    { // columns
      starts.read( reinterpret_cast<char*>(&startDoubles ), sizeof startDoubles );
      arrayCPU[j][i] = startDoubles;
      //std::cout<< arrayCPU[j][i]<<"\n";
    }
  }

  starts.close();

  double cost;

  for (int j = 0; j < size; j++)
  {
    for(int i = 0; i < parameterSize; i++)
    {
      singleArray[i]=arrayCPU[j][i]; // Dint knwo how to pass the double array as a 1D array to the trajectory Print function
    }
    double  lastStep;
    double  cost;
    int n = -1; // So it doesnt write files - if you want to rewrite the solution instead of n use i
    elements<double> yOut; // If you want to check the position and velocity of the last step
    cost = trajectoryPrint(singleArray,lastStep, cost, n, yOut);
    //writeTrajectoryToFile(singleArray, cost,j);
    std::cout<<std::endl<<"Run: "<< j+1 <<std::endl;
    std::cout<<cost<<std::endl;
    //std::cout<<yOut<<launchCon->getCondition(singleArray[TRIPTIME_OFFSET])<<std::endl;
  } 

}