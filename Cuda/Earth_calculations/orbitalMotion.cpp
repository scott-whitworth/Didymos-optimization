//Didymos-Optimization_Project:
//Last Editor: Mateo, Lauren, and Ben
//Tasks Completed: 
    //No recent changes

#include "../Runge_Kutta/runge_kutta.h" // used for rk4Reverse().
#include <iostream> // used for cout
#include <fstream> // stream output
#include <string> // to_string()
#include "earthInfo.h" // reference to launchCon

elements<double> earthInitial_incremental(double timeInitial, double tripTime,const elements<double> & earth) {
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

elements<double> earthInitial(double timeInitial, double tripTime,const elements<double> & earth) {
  // Time step
  double deltaT; 

  // Initial guess for time step, cannot be greater than the time resolution.
  deltaT = -(tripTime - timeInitial)/static_cast <double> (MAX_NUMSTEPS); 

  // Declaring the solution vector.
  elements<double> yp;

  // Calculates the earth's launch date conditions based on timeFinal minus the optimized trip time.
  rk4Reverse(timeInitial,tripTime,earth,deltaT,yp,RK_TOL);
 
  return yp;
}


//taken from CPU code to output final results of genetic algorithm
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------

double trajectoryPrint( double x[], double & lastStep, double & cost, int j, elements<double> & yOut, thruster<double> thrust, int runNumber) {
  /*set the asteroid and inital conditions for the earth and spacecraft:
  constructor takes in radial position(au), angluar position(rad), off-plane position(au),
  radial velocity(au/s), azimuthal velocity(rad/s), off-plane velocity(au/s)*/

  // setting landing conditions of the asteroid (October 5, 2022)
  elements<double> asteroid = elements<double>(R_FIN_AST, THETA_FIN_AST, Z_FIN_AST, VR_FIN_AST, VTHETA_FIN_AST, VZ_FIN_AST);

  // setting initial conditions of earth based off of the impact date (October 5, 2022) minus the trip time (optimized).
  elements<double> earth =  launchCon->getCondition(x[TRIPTIME_OFFSET]);
  
  // setting initial conditions of the spacecraft
  elements<double> spaceCraft = elements<double>(earth.r+ESOI*cos(x[ALPHA_OFFSET]), earth.theta + asin(sin(M_PI-x[ALPHA_OFFSET])*ESOI/earth.r), earth.z,
                                                 earth.vr + cos(x[ZETA_OFFSET])*sin(x[BETA_OFFSET])*vEscape, earth.vtheta + cos(x[ZETA_OFFSET])*cos(x[BETA_OFFSET])*vEscape,
                                                 earth.vz + sin(x[ZETA_OFFSET])*vEscape);

  // setting time parameters
  double timeInitial=0; 
  double timeFinal=orbitalPeriod; // Orbital period of asteroid(s)
  double deltaT; // time step
  int numSteps = 5000; // initial guess for the number of time steps, guess for the memory allocated 
  deltaT = (timeFinal-timeInitial) / MAX_NUMSTEPS; // initial guess for time step, small is preferable

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
  int lastStepInt;

  rk4sys(timeInitial, x[TRIPTIME_OFFSET] , times, spaceCraft, deltaT, yp, absTol, coeff, accel, gamma, tau, lastStepInt, accel_output, fuelSpent, wetMass, thrust);

  lastStep = lastStepInt;

  // gets the final y values of the spacecrafts for the cost function.
  yOut = yp[(int)lastStep];
  // cost equation determines how close a given run is to impact.
  // based off the position components of the spacecraft and asteroid.
  double cost_pos/*, vel*/;
  cost_pos = pow(asteroid.r-yOut.r,2)+pow(asteroid.theta-fmod(yOut.theta,2*M_PI),2)+pow(asteroid.z-yOut.z,2);         
//  vel = sqrt(pow(asteroid.vr-yOut.vr,2)+pow(asteroid.vtheta-yOut.vtheta,2)+pow(asteroid.vz-yOut.vz,2));
  //cost = cost_pos<cost_vel?cost_pos:cost_vel;
  cost = cost_pos;

  // Flag to not print the solution
  if (j > 0) {
    // when the cost function is less than 10^-20, it is set to 0 in order to keep that answer of optimized values.
    if (cost < Fmin) {
      cost = 0;
    }

    // Output of yp to a binary file
    std::ofstream output;
    
    output.open ("orbitalMotion-accel"+std::to_string(runNumber)+".bin", std::ios::binary); 

    for(int i = 0; i <= lastStep; i++) {
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

void writeTrajectoryToFile(double *start, double & cost, int i, thruster<double> thrust, int runNumber) {
  double numStep = 0; // must be double to match output to binary file
  elements<double> yp;
  trajectoryPrint(start, numStep, cost, i, yp, thrust, runNumber);
  //writes final optimization values to a seperate file
  std::ofstream output;

  output.open ("final-optimization"+std::to_string(runNumber)+".bin", std::ios::binary);

  for (int j = 0; j < OPTIM_VARS; j++) {
    output.write((char*)&start[j], sizeof (double));
  }

  output.write((char*)&numStep, sizeof (double));
  output.close();
}