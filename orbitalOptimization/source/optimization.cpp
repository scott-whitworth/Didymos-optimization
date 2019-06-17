#include <iostream> // cout
#include <iomanip> //used for setw(), sets spaces between values

#include "nelder_mead.h"
#include "orbitalMotion.h"

//  Licensing:
//    This code is distributed under the GNU LGPL license. 
//  Modified:
//    27 February 2008
//  Author:
//    John Burkardt


int main ();

// optimizes a vector for orbital motion and writes a binary file for the best solution 
void optimizing ();

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
  n = 14; 

  // allocating memory according to number of variables
  start = new double[n];
  step = new double[n];
  xmin = new double[n];

  std::cout << "\n"<<"minimizing orbital motion"<<std::endl;

  // x[0]-x[8]: gamma coefficients
  // x[9]-x[11]: tau coefficients
  // x[12]: alpha - launch angle (declination) position 
  // x[13]: beta - launch angle (declination) velocity 

  // initial guesses for variables
  start[0] = 1.5;
  start[1] = 1.5;
  start[2] = 1.5;
  start[3] = 1.5;
  start[4] = 1.5;
  start[5] = 1.5;
  start[6] = 1.5;
  start[7] = 1.5;
  start[8] = 1.5;
  start[9] = 1.5;
  start[10] = 1.5;
  start[11] = 1.5;
  start[12] = 0.5;
  start[13] = 0.5;

  // convergence tolerance
  reqmin = 1.0E-26;

  // initial change in variable size
  // based off of the variable start value
  step[0] = 1.0E02;
  step[1] = 1.0E02;
  step[2] = 1.0E02;
  step[3] = 1.0E02;
  step[4] = 1.0E02;
  step[5] = 1.0E02;
  step[6] = 1.0E02;
  step[7] = 1.0E02;
  step[8] = 1.0E02;
  step[9] = 1.0E02;
  step[10] = 1.0E02;
  step[11] = 1.0E02;
  step[12] = 1.0E00;
  step[13] = 1.0E00;

  // how often the equation checks for a convergence
  konvge = 15;
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
  std::cout << "\n"<< "  Return code IFAULT = " << ifault << "\n"<< "  Estimate of minimizing value X*:\n"<< "\n";
  for (i = 0; i < n; i++)
  {
    std::cout << std::setw(2) << xmin[i] << ",";
  }
  std::cout << "\n" << "  F(X*) = " << ynewlo << "\n";
  std::cout << "\n"<< "  Number of iterations = " << icount << "\n"<< "  Number of restarts =   " << numres << "\n";

  // writes the solution based on optimized variables to a binary file
  trajectoryPrint(xmin);

  // cleans up dynamic memory
  delete [] start;
  delete [] step;
  delete [] xmin;

  return;
}
