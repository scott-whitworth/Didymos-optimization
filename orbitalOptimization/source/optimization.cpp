#include "optimization.h" 
#include "nelder_mead.h" // used for nelmin()
#include "constants.h" //used for wetMass
#include "orbitalMotion.h" //used for trajectory() and trajectoryPrint()
#include <iostream> // cout
#include <iomanip> //used for setw(), sets spaces between values


//  Licensing:
//    This code is distributed under the GNU LGPL license. 
//  Modified:
//    27 February 2008
//  Author:
//    John Burkardt

int main ()
{
  timestamp ();
  std::cout << "\n"<<"beginning of optimization"<<std::endl;

  optimizing ();

  std::cout << "\n"<<"end of optimization"<<std::endl;
  timestamp ();
  
  return 0;
}

void optimizing ()
//  Purpose:
//    optimizes the coefficients for gamma and tau fourier series. Also, optimizes alpha and beta angles (used in initial velocity of the spacecraft).
{
  // initializing variables for nelmin algorithm. See nelder_mead.cpp for input/output information
  int i; 
  int icount; 
  int ifault;
  int kcount;
  int konvge;
  int n;
  int numres;
  double reqmin;
  double *start;
  double *step;
  double *xmin;
  double ynewlo;

  // number of variables to be optimized
  n = 22; 

  // allocating memory according to number of variables
  start = new double[n];
  step = new double[n];
  xmin = new double[n];

  std::cout << "\n"<<"minimizing orbital motion"<<std::endl;

  // x[0]-x[8]: gamma coefficients used to calculate fourier series
  // x[9]-x[11]: tau coefficients used to calculate fourier series
  // x[12]: alpha - launch angle (declination) position 
  // x[13]: beta - launch angle (declination) velocity 
  // x[14]: trip time - total time from launch to impact, sets the initial earth position
  // x[15-19]: coast coefficients used to calculate fourier series
  // x[20]: coast threshold - value set to determine when coasting occurs
  // x[21]: dry mass - mass of spacecraft excluding fuel

  // initial guesses for variables based off of previous runs which have small cost values
  start[0] = 10;
  start[1] = 10;
  start[2] = 10;
  start[3] = 10;
  start[4] = 10;
  start[5] = 10;
  start[6] = 10;
  start[7] = 10;
  start[8] = 10;
  start[9] = 10;
  start[10] = 10;
  start[11] = 10;
  start[12] = 0.5;
  start[13] = 0.5;
  start[14] = 6.653820100923719e+07/2; // period of asteroid divided by two, approx. 1.1 years
  start[15] = 1;
  start[16] = 1;
  start[17] = 1;
  start[18] = 1;
  start[19] = 1;
  start[20] = 0.25;
  start[21] = wetMass-250; // 2700 kg

  // convergence tolerance
  reqmin = 1.0E-40;

  // initial change in variable size
  // based off of the variable start value
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
  step[20] = 5.0E-01;
  step[21] = 5.0E01;

  // how often the equation checks for a convergence
  konvge = 20;
  // maximum number of iterations for convergence
  kcount = 10000;


  std::cout << "\n"<<"starting conditions"<<std::endl;
  for ( i = 0; i < n; i++ )
  {
    std::cout << std::setw(2) << start[i] << "\n";
  }

  // optimization value for the initial conditions
  ynewlo = trajectory (start);

  std::cout << "\n"<< " F(X) = " << ynewlo << std::endl;
  
  // nelder_mead function (optimization function)
  // see nelder_mead.cpp for input and output information
  nelmin (trajectory, n, start, xmin, &ynewlo, reqmin, step,
    konvge, kcount, &icount, &numres, &ifault);

  // displays error type when an error occurs
  std::cout << "\nReturn code IFAULT = " << ifault << "\nEstimate of minimizing value X*:\n\n";
  for (i = 0; i < n; i++)
  {
    std::cout << std::setw(2) << xmin[i] << ",";
  }
  std::cout << "\nF(X*) = " << ynewlo << "\n";
  std::cout << "\n"<< "  Number of iterations = " << icount << "\n"<< "  Number of restarts =   " << numres << "\n";

  // writes the solution based on optimized variables to a binary file
  trajectoryPrint(xmin);

  // cleans up dynamic memory
  delete [] start;
  delete [] step;
  delete [] xmin;

  return;
}
