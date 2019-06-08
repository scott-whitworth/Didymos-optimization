#include "rk4sys.h"
#include "calcThrust.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include <chrono>
#include <ratio>
//TODO: SC: Why are we using these? What benefit are they serving?

int main()
{
    // setting initial conditions of the asteroid
    elements<double> y0;
    y0.r = 3.150802646376772e+11/AU;// radial position (au)
    y0.theta= -3.081519548404041;// angular position (rad)
    y0.z =  1.760293325286572e+10/AU;// off-plane position (au)
    y0.vr = 4706.64912336045/AU;// radial velocity (au/s)
    y0.vtheta= 16716.9055348804/AU;// azimuthal velocity (rad/s)
    y0.vz= -81.4453413932308/AU;// off-plane velocity (au/s)


    double accel = 0.0001/AU;// thrust acceleration (au/s^2)


    // setting time parameters
    double timeInitial=0; 
    double timeFinal=6.653820100923719e+07; // Orbital period of asteroid(s)
    double deltaT; // time step
    int numSteps = 5000; // initial guess for the number of time steps
    deltaT = (timeFinal-timeInitial)/1e9; // initial guess for time step, small is preferable

    // setup of thrust angle calculations
    coefficients<double> coeff;
    for (int i=1;i<=10;i++){
      coeff.gamma[i]=0;
    }
    for (int i=1;i<=5;i++){
      coeff.tau[i]=0;
    }

    // setting Runge-Kutta tolerance
    double absTol = 1e-9;

    // Initialize memory for the solution vector of the dependant solution
    //TODO: SC: That is *times? Where and why are you using numSteps == 500?
    elements<double>* yp;
    yp = new elements<double>[numSteps];
    double *times;
    times = new double[numSteps];

  
    for (int repeat = 0; repeat<1; repeat++){
      yp = rk4sys(timeInitial,timeFinal,times,y0,deltaT,yp,absTol,coeff,accel);
    }
    std::cout<<times[1]<<yp[1]<<std::endl;

// Output of yp to a binary file
  std::ofstream output;
  
  output.open ("orbitalMotion-accel.bin", std::ios::binary);
  for(int i=0; i < numSteps; i++)
  {
    //output << yp[i];
    output.write((char*)&yp[i], sizeof (elements<double>));
    output.write((char*)&times[i], sizeof (double));
  }
  output.close();

  // cleaning up dynamic yp and time
    delete [] yp;
    delete [] times;

    return 0;
}