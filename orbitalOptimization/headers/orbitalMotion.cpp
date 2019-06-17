#include "rk4sys.h"
#include "calcThrust.h"
#include <iostream> // used for cout
#include <fstream> // used for stream output 
#include <ctime> // used for clock
#include <chrono> // used for clock
#include <math.h>

double trajectory( double x[])
{
// setting the acceleration as a constant (temporary)
    double accel = 0.001/AU;// thrust acceleration (au/s^2)

    // set landing conditions for Earth and the asteroid and inital conditions for the spacecraft:
    // constructor takes in radial position(au), angluar position(rad), off-plane position(au),
    // radial velocity(au/s), azimuthal velocity(rad/s), off-plane velocity(au/s)

    // setting landing conditions of the asteroid (October 5, 2022)
    elements<double> asteroid = elements<double>(1.02696822710421, 0.238839574416454, -0.0526614832914496,
    -2.05295246185041e-08, 2.29132593453064e-07, 8.00663905822009e-09);

    // setting landing conditions of earth (October 5, 2022)
    elements<double> earth = elements<double>(1.00021392223428, 0.199470650149394, -1.54878511585620e-05,
    -3.32034068725821e-09, 1.99029138292504e-07, -9.71518257891386e-12);

    // setting initial conditions of the spacecraft
    // not the actual initial conditions, right now just equal to the earth's landing date conditions
    elements<double> spaceCraft = elements<double>(earth.r+ESOI*cos(x[12]), earth.theta+asin(sin(M_PI-x[12])*ESOI/earth.r),earth.z,
    earth.vr+sin(x[13])*vEscape, earth.vtheta+cos(x[13])*vEscape,earth.vz);

    // setting time parameters
    double timeInitial=0; 
    double timeFinal=Torbital/5; // Orbital period of asteroid(s)
    double deltaT; // time step
    deltaT = (timeFinal-timeInitial)/1e9; // initial guess for time step, small is preferable

    // setup of thrust angle calculations
    coefficients<double> coeff;
    for (int i=0;i<coeff.gammaSize;i++){
      coeff.gamma[i]=x[i];
    }
    for (int i=0;i<coeff.tauSize;i++){
      coeff.tau[i]=x[i+9];
    }
    
    // setting Runge-Kutta tolerance
    double absTol = 1e-12;

    //set optmization minimum
    double Fmin = 1e-18;

    // Initialize memory for the solution vector of the dependant solution
    elements<double> yp;

    int lastStep = 0;
  
    rk4Simple(timeInitial,timeFinal,spaceCraft,deltaT,yp,absTol,coeff,accel,lastStep);

  double cost;
  cost = pow(asteroid.r-yp.r,2)+pow(asteroid.theta-yp.theta,2)+pow(asteroid.z-yp.z,2);

  if (cost < Fmin)
    cost = 0;
  
 std::cout<<"The cost value is: "<<cost<<std::endl;


  return cost;
}

double trajectoryPrint( double x[])
{
// setting the acceleration as a constant (temporary)
    double accel = 0.001/AU;// thrust acceleration (au/s^2)

    // set landing conditions for Earth and the asteroid and inital conditions for the spacecraft:
    // constructor takes in radial position(au), angluar position(rad), off-plane position(au),
    // radial velocity(au/s), azimuthal velocity(rad/s), off-plane velocity(au/s)

    // setting landing conditions of the asteroid (October 5, 2022)
    elements<double> asteroid = elements<double>(1.02696822710421, 0.238839574416454, -0.0526614832914496,
    -2.05295246185041e-08, 2.29132593453064e-07, 8.00663905822009e-09);

    // setting landing conditions of earth (October 5, 2022)
    elements<double> earth = elements<double>(1.00021392223428, 0.199470650149394, -1.54878511585620e-05,
    -3.32034068725821e-09, 1.99029138292504e-07, -9.71518257891386e-12);

    // setting initial conditions of the spacecraft
    // not the actual initial conditions, right now just equal to the earth's landing date conditions
    elements<double> spaceCraft = elements<double>(earth.r+ESOI*cos(x[12]), earth.theta+asin(sin(M_PI-x[12])*ESOI/earth.r),earth.z,
    earth.vr+sin(x[13])*vEscape, earth.vtheta+cos(x[13])*vEscape,earth.vz);

    // setting time parameters
    double timeInitial=0; 
    double timeFinal=Torbital/5; // Orbital period of asteroid(s)
    double deltaT; // time step
    int numSteps = 5000; // initial guess for the number of time steps, guess for the memory allocated 
    deltaT = (timeFinal-timeInitial)/1e9; // initial guess for time step, small is preferable

    // setup of thrust angle calculations
    coefficients<double> coeff;
    for (int i=0;i<coeff.gammaSize;i++){
      coeff.gamma[i]=x[i];
    }
    for (int i=0;i<coeff.tauSize;i++){
      coeff.tau[i]=x[i+9];
    }

    // setting Runge-Kutta tolerance
    double absTol = 1e-12;

    //set optmization minimum
    double Fmin = 1e-18;

    // Initialize memory for the solution vector of the dependant solution
    elements<double>* yp;
    yp = new elements<double>[numSteps];
    // Initialize memory for the values of the independent variable
    double *times;
    times = new double[numSteps];

    double *gamma;
    gamma = new double[numSteps];

    double *tau;
    tau = new double[numSteps];

    int lastStep = 0;
  
    rk4sys(timeInitial,timeFinal,times,spaceCraft,deltaT,yp,absTol,coeff,accel,gamma,tau,lastStep);

  elements<double> yFinal;
  yFinal = yp[lastStep];
  
  // cleaning up dynamic yp, time, gamma, and tau.


  double cost;
  cost = pow(asteroid.r-yFinal.r,2)+pow(asteroid.theta-yFinal.theta,2)+pow(asteroid.z-yFinal.z,2);

  if (cost < Fmin)
    cost = 0;
  
 std::cout<<"The cost value is: "<<cost<<std::endl<<yFinal;

  
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
  }
  output.close();
  

  delete [] yp;
  delete [] times;
  delete [] gamma;
  delete [] tau;

  return cost;
}

