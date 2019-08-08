// Didymos-Optimization_Project:
// Last Editor: Mateo and Lauren
// Tasks Completed: 
    // Added the zeta launch angle

#include "runge_kutta.h" // used for rk4sys(), rk4Simple90, and rk4Reverse().
#include "calcFourier.h" // used for calc_gamma(), calc_tau(), and calc_coast().
#include "earthInfo.h"
#include <iostream> // used for cout
#include <fstream> // used for stream output 
#include <math.h> // used for M_PI

// Solves orbital motion differential equations according to a vector of parameters (which are optimized) and returns the cost for the parameters
double trajectory( double x[])
{
  // Defining acceleration
    double accel;
    double massFuelSpent = 0;

  /*Set the asteroid and inital conditions for the earth and spacecraft:
  constructor takes in radial position(au), angluar position(rad), off-plane position(au),
  radial velocity(au/s), azimuthal velocity(rad/s), off-plane velocity(au/s)*/

  // Setting impact conditions of the asteroid (October 5, 2022)
  elements<double> asteroid = elements<double>(R_FIN_AST, THETA_FIN_AST, Z_FIN_AST, VR_FIN_AST, VTHETA_FIN_AST, VZ_FIN_AST);

  // Setting initial conditions of earth based on the impact date (October 5, 2022) minus the trip time (optimized)
  elements<double> earth = launchCon->getCondition(x[TRIPTIME_OFFSET]);
  // Setting initial conditions of the spacecraft
  elements<double> spaceCraft = elements<double>(earth.r+ESOI*cos(x[ALPHA_OFFSET]), //earth.r+ESOI*cos(alpha)

  /*Initial Angular Position is calculated using law of sines and the 
  triangle made by ESOI Radius, the radius of earth, and the radius from
  the Sun to the spacecraft. The angle at the suns apex is the
  calculated and added to the angle of the Earth to calculate the angle
  of the spacecraft.*/
  earth.theta+asin(sin(M_PI-x[ALPHA_OFFSET])*ESOI/earth.r), // earth.theta +arcsin(sin(pi-alpha))*ESOI/earth.r

  earth.z,

  /* Calculates initial radial velocity using earths radial velocity, the ship's
  scalar velocity, and angle of escape relative to Earth.
  escapeEarthAngle is the escape angle relative to earth and
  is defined to be 0 when the velocity is entirely angular
  and 90 when it is entirely radial. This equation calculates the ships
  initial radius from the sun by combining these values.*/
  earth.vr+cos(x[ZETA_OFFSET])*sin(x[BETA_OFFSET])*vEscape, // earth.vr + sin(beta)*vEscape

  /*Calculates initial specific angular momementum of ship using earth's
  specific angular momementum, the ships scalar velocity, escape angle,
  and initial radius.*/
  earth.vtheta+cos(x[ZETA_OFFSET])*cos(x[BETA_OFFSET])*vEscape, // earth.vtheta + cos(beta)*vEscape

  earth.vz+sin(x[ZETA_OFFSET])*vEscape);

  // setting time parameters
  double timeInitial=0; 
  double timeFinal=orbitalPeriod; // Orbital period of asteroid(s)
  double deltaT; // time step
  deltaT = (timeFinal-timeInitial)/MAX_NUMSTEPS; // initial guess for time step, small is preferable

  // Setup of thrust angle calculations based off of optimized coefficients
  // In-plane angle
  coefficients<double> coeff;
  initCoefficient(x,coeff);
  // Assigning coast threshold (now done in coefficients because is a constant)
 
  // Sssigning wetMass
  double wetMass = WET_MASS;

  // Setting Runge-Kutta tolerance
  double absTol = RK_TOL;

  // Set optmization minimum
  double Fmin = F_MIN;

  // Initialize memory for the solution vector of the dependant solution
  elements<double> yp;


  // Calling rk4simple for efficieny, calculates the trip data based on the final optimized value of y
  rk4Simple(timeInitial,x[TRIPTIME_OFFSET],spaceCraft,deltaT,yp,absTol,coeff,accel,wetMass);
  // Cost equation determines how close a given run is to impact.
  // Based off the position components of the spacecraft and asteroid.
  double cost, cost_pos, cost_vel;
  cost_pos = pow(asteroid.r-yp.r,2)+pow(asteroid.theta-yp.theta,2)+pow(asteroid.z-yp.z,2);
  cost_vel = pow((sqrt(pow(asteroid.vr-yp.vr,2)+pow(asteroid.vtheta-yp.vtheta,2)+pow(asteroid.vz-yp.vz,2))-V_IMPACT)/V_IMPACT,2);
  //cost = cost_pos<cost_vel?cost_pos:cost_vel;
  cost = cost_pos;

  //+pow((wetMass-DRY_MASS-massFuelSpent)/(wetMass-DRY_MASS),2)
  // when the cost function is less than 10^-20, it is set to 0 in order to keep that answer of optimized values.
  if (cost < Fmin)
    cost = 0;

  // output of the cost value
  //std::cout<<"The cost value is: "<<cost<<std::endl;

  return cost;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double trajectoryPrint( double x[], double & lastStep, double & cost, int j, elements<double> & yOut)
{
  /*set the asteroid and inital conditions for the earth and spacecraft:
  constructor takes in radial position(au), angluar position(rad), off-plane position(au),
  radial velocity(au/s), azimuthal velocity(rad/s), off-plane velocity(au/s)*/

  // setting landing conditions of the asteroid (October 5, 2022)
  elements<double> asteroid = elements<double>(R_FIN_AST, THETA_FIN_AST, Z_FIN_AST, VR_FIN_AST, VTHETA_FIN_AST, VZ_FIN_AST);

  // setting initial conditions of earth based off of the impact date (October 5, 2022) minus the trip time (optimized).
  elements<double> earth =  launchCon->getCondition(x[TRIPTIME_OFFSET]);
  
  // setting initial conditions of the spacecraft
  elements<double> spaceCraft = elements<double>(earth.r+ESOI*cos(x[ALPHA_OFFSET]), earth.theta+asin(sin(M_PI-x[ALPHA_OFFSET])*ESOI/earth.r),earth.z,
  earth.vr+cos(x[ZETA_OFFSET])*sin(x[BETA_OFFSET])*vEscape, earth.vtheta+cos(x[ZETA_OFFSET])*cos(x[BETA_OFFSET])*vEscape,earth.vz+sin(x[ZETA_OFFSET])*vEscape);

  // setting time parameters
  double timeInitial=0; 
  double timeFinal=orbitalPeriod; // Orbital period of asteroid(s)
  double deltaT; // time step
  int numSteps = 5000; // initial guess for the number of time steps, guess for the memory allocated 
  deltaT = (timeFinal-timeInitial)/MAX_NUMSTEPS; // initial guess for time step, small is preferable

  // setup of thrust angle calculations based off of optimized coefficients
  coefficients<double> coeff;
  initCoefficient(x,coeff);
  // Assigning coast threshold (now done in coefficients because is a constant)

  // Assigning wetMass
  double wetMass = WET_MASS;
  // setting Runge-Kutta tolerance
  double absTol = RK_TOL;
  //set optmization minimum
  double Fmin = F_MIN;

  // Initialize memory for the solution vector of the dependant solution
  elements<double>* yp;
  yp = new elements<double>[numSteps];
  
  double *times, *gamma, *tau, *accel_output, *fuelSpent;
  times = new double[numSteps]; // Initialize memory for time array
  gamma = new double[numSteps]; // Initialize memory for gamma array
  tau = new double[numSteps]; // Initialize memory for tau array
  accel_output = new double[numSteps]; // Initialize memory for acceleration array
  fuelSpent = new double[numSteps];  // Initialize memory for fuelSpent array

  double accel; // Initialize memory for  acceleration

  // used to track the cost function throughout a run via output and outputs to a binary
  rk4sys(timeInitial,x[TRIPTIME_OFFSET],times,spaceCraft,deltaT,yp,absTol,coeff,accel,gamma,tau,lastStep,accel_output,fuelSpent, wetMass);

  // gets the final y values of the spacecrafts for the cost function.
  yOut = yp[(int)lastStep];
  // cost equation determines how close a given run is to impact.
  // based off the position components of the spacecraft and asteroid.
  double cost_pos, vel;
  cost_pos = pow(asteroid.r-yOut.r,2)+pow(asteroid.theta-fmod(yOut.theta,2*M_PI),2)+pow(asteroid.z-yOut.z,2);         
  vel = sqrt(pow(asteroid.vr-yOut.vr,2)+pow(asteroid.vtheta-yOut.vtheta,2)+pow(asteroid.vz-yOut.vz,2));
  //cost = cost_pos<cost_vel?cost_pos:cost_vel;
  cost = cost_pos;

  // Flag to not print the solution
  if (j>0)
  {
    // when the cost function is less than 10^-20, it is set to 0 in order to keep that answer of optimized values.
    if (cost < Fmin)
      cost = 0;

    std::cout<<"Impact velocity: "<<vel*AU<<" m/s"<<std::endl;
    // Output of yp to a binary file
    std::ofstream output;
    
    output.open ("orbitalMotion-accel"+std::to_string(j)+".bin", std::ios::binary); 

    for(int i=0; i <= lastStep; i++)
    {
      //output << yp[i];
      output.write((char*)&yp[i], sizeof (elements<double>));
      output.write((char*)&times[i], sizeof (double));
      output.write((char*)&gamma[i], sizeof (double));
      output.write((char*)&tau[i], sizeof (double));
      output.write((char*)&accel_output[i], sizeof (double));
      output.write((char*)&fuelSpent[i], sizeof (double));
    }
    output.close();
  }
  // cleaning up dynamic yp, time, gamma, and tau.
  delete [] yp;
  delete [] times;
  delete [] gamma;
  delete [] tau;

  return cost;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Solver the integration backwards but it is used in earth info, it is meant to be called consecutively
elements<double> earthInitial_incremental(double timeInitial, double tripTime,const elements<double> & earth)
{
  // Time step
  double deltaT; 

  // Initial guess for time step, cannot be greater than the time resolution.
  deltaT = -(tripTime - timeInitial)/static_cast <double> (60); 

  // Declaring the solution vector.
  elements<double> yp;

  // Calculates the earth's launch date conditions based on timeFinal minus the optimized trip time.
  rk4Reverse(timeInitial,tripTime,earth,deltaT,yp,RK_TOL);
 
  return yp;
}

// solves orbital motion differential equations according to a vector of parameters (which are optimized) and returns the cost for the parameters
// reverse integration in order to determine the initial conditions of the earth (at launch)
elements<double> earthInitial(double tripTime)
{
  //setting initial conditions for calculation of earth on launch date with orbital elements of the earth on the asteroid impact date of 2022-10-05.
  elements<double> earth = elements<double>(R_FIN_EARTH, THETA_FIN_EARTH, Z_FIN_EARTH, VR_FIN_EARTH, VTHETA_FIN_EARTH, VZ_FIN_EARTH);

  // setting intiial time parameters
  double timeInitial= 0; 
  double deltaT; // time step
  deltaT = -(tripTime-timeInitial)/static_cast <double> (MAX_NUMSTEPS); // initial guess for time step, small is preferable

  // declaring the solution vector
  elements<double> yp;

  // setting Runge-Kutta tolerance
  double absTol = RK_TOL;

  // calculates the earth's launch date conditions based on timeFinal minus the optimized trip time
  rk4Reverse(timeInitial,tripTime,earth,deltaT,yp,absTol);

  return yp;
}


