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

  n = 14;

  start = new double[n];
  step = new double[n];
  xmin = new double[n];

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Apply NELMIN to ROSENBROCK function.\n";

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

  reqmin = 1.0E-26;

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
  step[12] = 1.0E01;
  step[13] = 1.0E01;

  konvge = 15;
  kcount = 10000;

  cout << "\n";
  cout << "  Starting point X:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "" << setw(14) << start[i] << "\n";
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
    cout << setw(14) << xmin[i] << ",";
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