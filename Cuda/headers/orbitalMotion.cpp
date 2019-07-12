//Didymos-Optimization_Project:
//Last Editor: Mateo, Lauren, and Ben
//Tasks Completed: 
    //No recent changes
#define _USE_MATH_DEFINES // for use of M_PI

#include "runge_kutta.h" // used for rk4sys(), rk4Simple(), and rk4Reverse().
#include "runge_kuttaCUDA.cuh" // used for rk4SimpleCUDA()
#include "calcFourier.h" // used for calc_gamma(), calc_tau(), and calc_coast().
#include <iostream> // used for cout
#include <fstream> // used for stream output 
#include <math.h> // used for M_PI

// Solves orbital motion differential equations according to a vector of parameters (which are optimized) and returns the cost for the parameters
double trajectory( double x[])
{
  // Defining acceleration
    double accel;
    //double massFuelSpent = 0;

  /*Set the asteroid and inital conditions for the earth and spacecraft:
  constructor takes in radial position(au), angluar position(rad), off-plane position(au),
  radial velocity(au/s), azimuthal velocity(rad/s), off-plane velocity(au/s)*/

  // Setting impact conditions of the asteroid (October 5, 2022)
  elements<double> asteroid = elements<double>(R_FIN_AST, THETA_FIN_AST, Z_FIN_AST, VR_FIN_AST, VTHETA_FIN_AST, VZ_FIN_AST);

  // Setting initial conditions of earth based on the impact date (October 5, 2022) minus the trip time (optimized)
  elements<double> earth =  earthInitial(x[TRIPTIME_OFFSET]);
  
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
  earth.vr+sin(x[BETA_OFFSET])*vEscape, // earth.vr + sin(beta)*vEscape

  /*Calculates initial specific angular momementum of ship using earth's
  specific angular momementum, the ships scalar velocity, escape angle,
  and initial radius.*/
  earth.vtheta+cos(x[BETA_OFFSET])*vEscape, // earth.vtheta + cos(beta)*vEscape

  earth.vz);

  // setting time parameters
  double timeInitial=0; 
  double timeFinal=orbitalPeriod; // Orbital period of asteroid(s)
  double deltaT; // time step
  deltaT = (timeFinal-timeInitial)/MAX_NUMSTEPS; // initial guess for time step, small is preferable

  // Setup of thrust angle calculations based off of optimized coefficients
  // In-plane angle
  coefficients<double> coeff;
  for (int i=0;i<coeff.gammaSize;i++){
    coeff.gamma[i]=x[i+GAMMA_OFFSET];
  }
  // Out-of-plane angle
  for (int i=0;i<coeff.tauSize;i++){
    coeff.tau[i]=x[i+TAU_OFFSET];
  }

  // Setup of coast determination calculations based off of optimized coefficients
  for (int i=0;i<coeff.coastSize;i++){
    coeff.coast[i]=x[i+COAST_OFFSET];
  //coeff.coast[i]=0.;

  }
  // Assigning optimized coast threshold
  coeff.coastThreshold = x[THRESHOLD_OFFSET];
  // coeff.coastThreshold = -1.;
 
  // Sssigning optimized wetMass
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
  double cost, cost_pos;
  //double cost_vel;
  cost_pos = pow(asteroid.r-yp.r,2)+pow(asteroid.theta-yp.theta,2)+pow(asteroid.z-yp.z,2);
  //cost_vel = pow((sqrt(pow(asteroid.vr-yp.vr,2)+pow(asteroid.vtheta-yp.vtheta,2)+pow(asteroid.vz-yp.vz,2))-V_IMPACT)/V_IMPACT,2);
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

double trajectoryPrint( double x[], int & n, double & cost)
{
  // defining the acceleration
    double accel;

  /*set the asteroid and inital conditions for the earth and spacecraft:
  constructor takes in radial position(au), angluar position(rad), off-plane position(au),
  radial velocity(au/s), azimuthal velocity(rad/s), off-plane velocity(au/s)*/

  // setting landing conditions of the asteroid (October 5, 2022)
  elements<double> asteroid = elements<double>(R_FIN_AST, THETA_FIN_AST, Z_FIN_AST, VR_FIN_AST, VTHETA_FIN_AST, VZ_FIN_AST);

  // setting initial conditions of earth based off of the impact date (October 5, 2022) minus the trip time (optimized).
  elements<double> earth =  earthInitial(x[TRIPTIME_OFFSET]);

  // setting initial conditions of the spacecraft
  elements<double> spaceCraft = elements<double>(earth.r+ESOI*cos(x[ALPHA_OFFSET]), earth.theta+asin(sin(M_PI-x[ALPHA_OFFSET])*ESOI/earth.r),earth.z,
  earth.vr+sin(x[BETA_OFFSET])*vEscape, earth.vtheta+cos(x[BETA_OFFSET])*vEscape,earth.vz);

  // setting time parameters
  double timeInitial=0; 
  double timeFinal=orbitalPeriod; // Orbital period of asteroid(s)
  double deltaT; // time step
  int numSteps = 5000; // initial guess for the number of time steps, guess for the memory allocated 
  deltaT = (timeFinal-timeInitial)/MAX_NUMSTEPS; // initial guess for time step, small is preferable

  // setup of thrust angle calculations based off of optimized coefficients
  // in-plane angle
  coefficients<double> coeff;
  for (int i=0;i<coeff.gammaSize;i++){
    coeff.gamma[i]=x[i+GAMMA_OFFSET];
  }
  // out-of-plane angle
  for (int i=0;i<coeff.tauSize;i++){
    coeff.tau[i]=x[i+TAU_OFFSET];
  }

  // setup of coast determination calculations based off of optimized coefficients
  for (int i=0;i<coeff.coastSize;i++){
    coeff.coast[i]=x[i+COAST_OFFSET];
  //  coeff.coast[i]=0.;
  }
  // assigning optimized coast threshold
  coeff.coastThreshold = x[THRESHOLD_OFFSET];
  //coeff.coastThreshold = -1.;
 
  // assigning optimized wetMass
  double wetMass = WET_MASS;
  // setting a resonable range for wetMass
  if(wetMass<DRY_MASS|| wetMass>3000)
  {
      return 100;
  }

  // setting Runge-Kutta tolerance
  double absTol = RK_TOL;

  //set optmization minimum
  double Fmin = F_MIN;

  // Initialize memory for the solution vector of the dependant solution
  elements<double>* yp;
  yp = new elements<double>[numSteps];
  // Initialize memory for time array
  double *times;
  times = new double[numSteps];
  // Initialize memory for gamma array
  double *gamma;
  gamma = new double[numSteps];
  // Initialize memory for tau array
  double *tau;
  tau = new double[numSteps];
  // Initialize memory for acceleration array
  double *accel_output;
  accel_output = new double[numSteps];

  // used to get yFinal
  int lastStep = 0;

  // used to track the cost function throughout a run via output and outputs to a binary
  rk4sys(timeInitial,x[TRIPTIME_OFFSET],times,spaceCraft,deltaT,yp,absTol,coeff,accel,gamma,tau,lastStep,accel_output, wetMass);

  // gets the final y values of the spacecrafts for the cost function.
  elements<double> yFinal;
  yFinal = yp[lastStep];
 
  // cost equation determines how close a given run is to impact.
  // based off the position components of the spacecraft and asteroid.
  double cost_pos;
  //double cost_vel;
  cost_pos = pow(asteroid.r-yFinal.r,2)+pow(asteroid.theta-yFinal.theta,2)+pow(asteroid.z-yFinal.z,2);         
  //cost_vel = pow((sqrt(pow(asteroid.vr-yFinal.vr,2)+pow(asteroid.vtheta-yFinal.vtheta,2)+pow(asteroid.vz-yFinal.vz,2))-V_IMPACT)/V_IMPACT,2);
  //cost = cost_pos<cost_vel?cost_pos:cost_vel;
  cost = cost_pos;

  // when the cost function is less than 10^-20, it is set to 0 in order to keep that answer of optimized values.
   if (cost < Fmin)
    cost = 0;

  n = lastStep;

  // Output of yp to a binary file
  std::ofstream output;
  
  output.open ("orbitalMotion-accel.bin", std::ios::binary); 

  for(int i=0; i <= lastStep; i++)
  {
    //output << yp[i];
    output.write((char*)&yp[i], sizeof (elements<double>));
    output.write((char*)&times[i], sizeof (double));
    output.write((char*)&gamma[i], sizeof (double));
    output.write((char*)&tau[i], sizeof (double));
    output.write((char*)&accel_output[i], sizeof (double));
  }
  output.close();

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

// solves orbital motion differential equations according to a vector of parameters (which are optimized) and returns the cost for the parameters
// reverse integration in order to determine the initial conditions of the earth (at launch)
elements<double> earthInitial(double tripTime)
{
  // setup of thrust angle calculations based off of optimized coefficients
  // all set to zero because the earth has no acceleration due to thrusting.
  coefficients<double> earthCoeff;
  for (int i=0;i<earthCoeff.gammaSize;i++){
    earthCoeff.gamma[i]=0;
  }
  for (int i=0;i<earthCoeff.tauSize;i++){
    earthCoeff.tau[i]=0;
  }

  // set to zero because the earth has no acceleration due to thrusting.
  double earthAccel = 0;

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
  rk4Reverse(timeInitial,tripTime,earth,deltaT,yp,absTol,earthCoeff,earthAccel);

  return yp;
}