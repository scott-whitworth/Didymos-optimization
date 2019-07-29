//Didymos-Optimization_Project:
//Last Editor: Mateo, Lauren, and Ben
//Tasks Completed: 
    //No recent changes

#include "../Runge_Kutta/runge_kutta.h" // used for rk4Reverse().
#include <iostream> // used for cout

elements<double> earthInitial_incremental(double timeInitial, double tripTime,const elements<double> & earth)
{

  //setting initial conditions for calculation of earth on launch date with orbital elements of the earth on the asteroid impact date of 2022-10-05.
  // setting intiial time parameters
  double deltaT; // time step
  deltaT = -(tripTime - timeInitial)/60; // initial guess for time step, small is preferable
  // declaring the solution vector
  elements<double> yp;
  // setting Runge-Kutta tolerance
  double absTol = RK_TOL;

  // calculates the earth's launch date conditions based on timeFinal minus the optimized trip time
  rk4Reverse(timeInitial,tripTime,earth,deltaT,yp,absTol);
 
  return yp;
}