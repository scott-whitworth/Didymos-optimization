#include "runge_kutta.h" 

#include <iostream> // used for cout
#include <fstream> // used for stream output 
#include <math.h> // used for M_PI


// solves orbital motion differential equations according to a vector of parameters (which are optimized) and returns the cost for the parameters
// reverse integration in order to determine the initial conditions of the earth (at launch)
elements<double> earthInitial(double timeInitial, double tripTime,const elements<double> & earth)
{

  //setting initial conditions for calculation of earth on launch date with orbital elements of the earth on the asteroid impact date of 2022-10-05.
 

  // setting intiial time parameters

  double deltaT; // time step
  deltaT = -(tripTime-timeInitial)/60; // initial guess for time step, small is preferable
  //deltaT = -3600;

  // declaring the solution vector
  elements<double> yp;

  // setting Runge-Kutta tolerance
  double absTol = RK_TOL;

  // calculates the earth's launch date conditions based on timeFinal minus the optimized trip time
  rk4Reverse(timeInitial,tripTime,earth,deltaT,yp,absTol);
 
  return yp;
}