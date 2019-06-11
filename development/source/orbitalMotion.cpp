#include "rk4sys.h"
#include "calcThrust.h"
#include <iostream> // used for cout
#include <fstream> // used for stream output 
#include <ctime> // used for clock
#include <chrono> // used for clock

int main()
{
    // setting landing conditions of the asteroid (October 5, 2022)
    elements<double> asteroid;
    asteroid.r = 1.02696822710421;// radial position (au)
    asteroid.theta= 0.238839574416454;// angular position (rad)
    asteroid.z = -0.0526614832914496;// off-plane position (au)
    asteroid.vr = -2.05295246185041e-08;// radial velocity (au/s)
    asteroid.vtheta= 2.29132593453064e-07;// azimuthal velocity (rad/s)
    asteroid.vz= 8.00663905822009e-09;// off-plane velocity (au/s)

    // setting landing conditions of earth (October 5, 2022)
    elements<double> earth;
    earth.r = 1.00021392223428;// radial position (au)
    earth.theta= 0.199470650149394;// angular position (rad)
    earth.z =  -1.54878511585620e-05;// off-plane position (au)
    earth.vr = -3.32034068725821e-09;// radial velocity (au/s)
    earth.vtheta= 1.99029138292504e-07;// azimuthal velocity (rad/s)
    earth.vz= -9.71518257891386e-12;// off-plane velocity (au/s)

    // setting initial conditions of the spacecraft
    // not the actual initial conditions, right now just equal to the earth's landing date conditions
    elements<double> spaceCraft;
    spaceCraft.r = earth.r;// radial position (au)
    spaceCraft.theta= earth.theta;// angular position (rad)
    spaceCraft.z = earth.z;// off-plane position (au)
    spaceCraft.vr = earth.vr;// radial velocity (au/s)
    spaceCraft.vtheta= earth.vtheta;// azimuthal velocity (rad/s)
    spaceCraft.vz=earth.vz;// off-plane velocity (au/s)

    // setting the acceleration as a constant (temporary)
    double accel = 0.0000/AU;// thrust acceleration (au/s^2)


    // setting time parameters
    double timeInitial=0; 
    double timeFinal=6.653820100923719e+07; // Orbital period of asteroid(s)
    double deltaT; // time step
    int numSteps = 5000; // initial guess for the number of time steps, guess for the memory allocated 
    deltaT = (timeFinal-timeInitial)/1e9; // initial guess for time step, small is preferable

    // setup of thrust angle calculations
    coefficients<double> coeff;
    for (int i=0;i<coeff.gammaSize;i++){
      coeff.gamma[i]=1;
    }
    for (int i=0;i<coeff.tauSize;i++){
      coeff.tau[i]=1;
    }
    // setting Runge-Kutta tolerance
    double absTol = 1e-12;

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

  
    // Recording the start time for performance metric
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    for (int repeat = 0; repeat<1; repeat++){
      yp = rk4sys(timeInitial,timeFinal,times,spaceCraft,deltaT,yp,absTol,coeff,accel,gamma,tau);
    }
     // recording stop time
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    // calculating elapsed time of rk4sys() call(s)
    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    std::cout << "rk4sys() call took " << time_span.count() << " seconds." << std::endl;

// Output of yp to a binary file
  std::ofstream output;
  
  output.open ("orbitalMotion-accel.bin", std::ios::binary);
  for(int i=0; i < numSteps; i++)
  {
    //output << yp[i];
    output.write((char*)&yp[i], sizeof (elements<double>));
    output.write((char*)&times[i], sizeof (double));
    output.write((char*)&gamma[i], sizeof (double));
    output.write((char*)&tau[i], sizeof (double));
  }
  output.close();


  // cleaning up dynamic yp, time, gamma, and tau.
    delete [] yp;
    delete [] times;
    delete [] gamma;
    delete [] tau;

    return 0;
}