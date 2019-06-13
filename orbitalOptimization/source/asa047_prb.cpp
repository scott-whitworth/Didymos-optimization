# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

#include "nelder_mead.h"
#include "orbitalMotion.h"

using namespace std;

int main ( );
void test01 ( );
double rosenbrock ( double x[2] );
void test04 ( );
double quartic ( double x[10] );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//    MAIN is the main program for ASA047_PRB.
//  Discussion:
//    ASA047_PRB tests the ASA047 library.
//  Licensing:
//    This code is distributed under the GNU LGPL license. 
//  Modified:
//    27 February 2008
//  Author:
//    John Burkardt
//
{
  double x[3];
 
  x[0] = -3.32034068725821e-09;
  x[1] = 1.99029138292504e-07;
  x[2] = -9.71518257891386e-12;
  trajectory(x);
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
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//    TEST01 demonstrates the use of NELMIN on 
//  Licensing:
//    This code is distributed under the GNU LGPL license. 
//  Modified:
//    27 February 2008
//  Author:
//    John Burkardt
//
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

  reqmin = 1.0E-9;

  step[0] = 1.0e-9;
  step[1] = 1.0e-10;
  step[2] = 1.0e-9;

  konvge = 10;
  kcount = 10000;

  cout << "\n";
  cout << "  Starting point X:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(14) << start[i] << "\n";
  }

  ynewlo = trajectoryP ( start );

  cout << "\n";
  cout << "  F(X) = " << ynewlo << "\n";

  nelmin ( trajectoryP, n, start, xmin, &ynewlo, reqmin, step,
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

  trajectoryP(xmin);

  delete [] start;
  delete [] step;
  delete [] xmin;

  return;
}
//****************************************************************************80

double rosenbrock ( double x[2] )

//****************************************************************************80
//
//  Purpose:
//
//    ROSENBROCK evaluates the Rosenbrock parabolic value function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 February 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    R ONeill,
//    Algorithm AS 47:
//    Function Minimization Using a Simplex Procedure,
//    Applied Statistics,
//    Volume 20, Number 3, 1971, pages 338-345.
//
//  Parameters:
//
//    Input, double X[2], the argument.
//
//    Output, double ROSENBROCK, the value of the function.
//
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
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//    TEST04 demonstrates the use of NELMIN on QUARTIC.
//  Licensing:
//    This code is distributed under the GNU LGPL license. 
//  Modified:
//    19 February 2008
//  Author:
//    John Burkardt
//
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
//****************************************************************************80

double quartic ( double x[10] )

//****************************************************************************80
//
//  Purpose:
//
//    QUARTIC evaluates a function defined by a sum of fourth powers.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 February 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    R ONeill,
//    Algorithm AS 47:
//    Function Minimization Using a Simplex Procedure,
//    Applied Statistics,
//    Volume 20, Number 3, 1971, pages 338-345.
//
//  Parameters:
//
//    Input, double X[10], the argument.
//
//    Output, double QUARTIC, the value of the function.
//
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
