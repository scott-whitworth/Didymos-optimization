#include "rk4sys.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include <chrono>
#include <ratio>


int main()
{
// setting initial conditions of the asteroid
elements<long double> y0;
y0.r = 3.150802646376772e+11/AU;// radial position (au)
y0.theta= -3.081519548404041;// angular position (rad)
y0.z =  1.760293325286572e+10/AU;// off-plane position (au)
y0.vr = 4706.64912336045/AU;// radial velocity (au/s)
y0.vtheta= 16716.9055348804/AU;// azimuthal velocity (rad/s)
y0.vz= -81.4453413932308/AU;// off-plane velocity (au/s)

// setting time parameters
long double timeInitial=0; 
long double timeFinal=6.653820100923719e+07; // Orbital period of asteroid(s)
long double deltaT; // time step
int numSteps = 500; // number of iterations
deltaT = (timeFinal-timeInitial)/numSteps;

// Initialize memory

// Initialize memory for the solution vector of the dependant solution
elements<long double>* y;
y = new elements<long double>[numSteps];

//TODO: verify and/or increase resolution of timing 

elements<long double> *yp;

// Comparing completion time to MATLAB's time
//    Recording the start time
std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

//    Our output function (yp)
yp = rk4sys<long double>(timeInitial,timeFinal,y0,deltaT,numSteps,y);
  
//    recording stop time
std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

//    calculating elapsed time
std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

//    Output of elasped time
std::cout << "It took me " << time_span.count() << "seconds." << std::endl;
 
// Output of yp to a text file
std::ofstream output;
  
output.open ("orbitalMotionLong.bin", std::ios::binary);
for(int i=0; i < numSteps; i++)
{
  output.write((char*)&yp[i], sizeof (elements<long double>));
}
output.close();

// cleaning up dynamic yp
delete [] yp;

return 0;
}