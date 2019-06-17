# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

#include "nelder_mead.h"
#include "orbitalMotion.h"

//  Licensing:
//    This code is distributed under the GNU LGPL license. 
//  Modified:
//    27 February 2008
//  Author:
//    John Burkardt
//****************************************************************************

using namespace std;

int main ( );
void test01 ( );
double rosenbrock ( double x[2] );
void test04 ( );
double quartic ( double x[10] );

int main ( )
{
  double x[11];
 
  timestamp ( );
  cout << "\n";
  cout << "ASA047_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the ASA047 library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "ASA047_PRB:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************

void test01 ( )

//****************************************************************************
//  Purpose:
//    TEST01 demonstrates the use of NELMIN on 

{
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

  n = 3;

  start = new double[n];
  step = new double[n];
  xmin = new double[n];

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Apply NELMIN to ROSENBROCK function.\n";

  start[0] = -3.32034068725821e-09;
  start[1] = 1.99029138292504e-07;
  start[2] = -9.71518257891386e-12;


  reqmin = 1.0E-10;

  step[0] = 1.0e-9;
  step[1] = 1.0e-10;
  step[2] = 1.0e-10;


  konvge = 15;
  kcount = 30000;

  cout << "\n";
  cout << "  Starting point X:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(14) << start[i] << "\n";
  }

  ynewlo = trajectory ( start );

  cout << "\n";
  cout << "  F(X) = " << ynewlo << "\n";

  nelmin ( trajectory, n, start, xmin, &ynewlo, reqmin, step,
    konvge, kcount, &icount, &numres, &ifault );

  cout << "\n";
  cout << "  Return code IFAULT = " << ifault << "\n";
  cout << "\n";
  cout << "  Estimate of minimizing value X*:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(14) << xmin[i] << "\n";
  }

  cout << "\n";
  cout << "  F(X*) = " << ynewlo << "\n";

  cout << "\n";
  cout << "  Number of iterations = " << icount << "\n";
  cout << "  Number of restarts =   " << numres << "\n";

  // One more time
  trajectoryPrint(xmin);

  delete [] start;
  delete [] step;
  delete [] xmin;

  return;
}
//****************************************************************************

double rosenbrock ( double x[2] )

//****************************************************************************
//  Purpose:
//    ROSENBROCK evaluates the Rosenbrock parabolic value function.
//  Parameters
//    Input, double X[2], the argument.
//    Output, double ROSENBROCK, the value of the function.

{
  double fx;
  double fx1;
  double fx2;

  fx1 = x[1] - x[0] * x[0];
  fx2 = 1.0 - x[0];

  fx = 100.0 * fx1 * fx1
     +         fx2 * fx2;

  return fx;
}
//****************************************************************************

void test04 ( )

//****************************************************************************
//  Purpose:
//    TEST04 demonstrates the use of NELMIN on QUARTIC.

{
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

  n = 10;

  start = new double[n];
  step = new double[n];
  xmin = new double[n];

  cout << "\n";
  cout << "TEST04\n";
  cout << "  Apply NELMIN to the QUARTIC function.\n";

  for ( i = 0; i < n; i++ )
  {
    start[i] = 1.0;
  }

  reqmin = 1.0E-08;

  for ( i = 0; i < n; i++ )
  {
    step[i] = 1.0;
  }

  konvge = 10;
  kcount = 500;

  cout << "\n";
  cout << "  Starting point X:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(14) << start[i] << "\n";
  }

  ynewlo = quartic ( start );

  cout << "\n";
  cout << "  F(X) = " << ynewlo << "\n";

  nelmin ( quartic, n, start, xmin, &ynewlo, reqmin, step,
    konvge, kcount, &icount, &numres, &ifault );

  cout << "\n";
  cout << "  Return code IFAULT = " << ifault << "\n";
  cout << "\n";
  cout << "  Estimate of minimizing value X*:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(14) << xmin[i] << "\n";
  }

  cout << "\n";
  cout << "  F(X*) = " << ynewlo << "\n";

  cout << "\n";
  cout << "  Number of iterations = " << icount << "\n";
  cout << "  Number of restarts =   " << numres << "\n";

  delete [] start;
  delete [] step;
  delete [] xmin;

  return;
}
//****************************************************************************

double quartic ( double x[10] )

//****************************************************************************
//  Purpose:
//    QUARTIC evaluates a function defined by a sum of fourth powers.
//  Parameters:
//    Input, double X[10], the argument.
//    Output, double QUARTIC, the value of the function.

{
  double fx;
  int i;

  fx = 0.0;

  for ( i = 0; i < 10; i++ )
  {
    fx = fx + x[i] * x[i] * x[i] * x[i];
  }

  return fx;
}
