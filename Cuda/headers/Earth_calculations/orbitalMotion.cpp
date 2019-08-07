//Didymos-Optimization_Project:
//Last Editor: Mateo, Lauren, and Ben
//Tasks Completed: 
    //No recent changes

#include "../Runge_Kutta/runge_kutta.h" // used for rk4Reverse().
#include <iostream> // used for cout

elements<double> earthInitial_incremental(double timeInitial, double tripTime,const elements<double> & earth)
{
  // Time step
  double deltaT; 

  // Initial guess for time step, cannot be greater than the time resolution.
  deltaT = -(tripTime - timeInitial)/60; 

  // Declaring the solution vector.
  elements<double> yp;

  // Calculates the earth's launch date conditions based on timeFinal minus the optimized trip time.
  rk4Reverse(timeInitial,tripTime,earth,deltaT,yp,RK_TOL);
 
  return yp;
}

elements<double> earthInitial(double timeInitial, double tripTime,const elements<double> & earth)
{
  // Time step
  double deltaT; 

  // Initial guess for time step, cannot be greater than the time resolution.
  deltaT = -(tripTime - timeInitial)/MAX_NUMSTEPS; 

  // Declaring the solution vector.
  elements<double> yp;

  // Calculates the earth's launch date conditions based on timeFinal minus the optimized trip time.
  rk4Reverse(timeInitial,tripTime,earth,deltaT,yp,RK_TOL);
 
  return yp;
}