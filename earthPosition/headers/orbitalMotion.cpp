//Didymos-Optimization_Project:
//Last Editor: Mateo, Lauren, and Ben
//Tasks Completed: 
    //No recent changes

#include "runge_kutta.h" // used for rk4sys(), rk4Simple90, and rk4Reverse().
#include <iostream> // used for cout
#include <fstream> // used for stream output 
#include <math.h> // used for M_PI


// solves orbital motion differential equations according to a vector of parameters (which are optimized) and returns the cost for the parameters
// reverse integration in order to determine the initial conditions of the earth (at launch)
elements<double> earthInitial(double tripTime)
{

  //setting initial conditions for calculation of earth on launch date with orbital elements of the earth on the asteroid impact date of 2022-10-05.
  elements<double> earth = elements<double>(R_FIN_EARTH, THETA_FIN_EARTH, Z_FIN_EARTH, VR_FIN_EARTH, VTHETA_FIN_EARTH, VZ_FIN_EARTH);

  // setting intiial time parameters
  double timeInitial= 0; 
  double deltaT; // time step
  deltaT = -(tripTime-timeInitial)/MAX_NUMSTEPS; // initial guess for time step, small is preferable

  // declaring the solution vector
  elements<double> yp;

  // setting Runge-Kutta tolerance
  double absTol = RK_TOL;

  // calculates the earth's launch date conditions based on timeFinal minus the optimized trip time
  rk4Reverse(timeInitial,tripTime,earth,deltaT,yp,absTol);

  return yp;
}