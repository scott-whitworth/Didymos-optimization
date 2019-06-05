#include "rk4sys.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include <chrono>
#include <ratio>


int main()
{
elements y0;
y0.r = 3.150802646376772e+11/AU;//Initial radial position of the asteroid (m)
y0.theta= -3.081519548404041;//Initial angular position of the asteroid (rad)
y0.z =  1.760293325286572e+10/AU;//Initial off-plane position of the asteroid (m)
y0.vr = 4706.64912336045/AU;//Initial radial velocity of the asteroid (m/s)
y0.vtheta= 16716.9055348804/AU;//Intial azimuthal velocity of the asteroid (m/s)
y0.vz= -81.4453413932308/AU;//Initial off-plane velocity of the asteroid (m/s)
double timeInitial=0; 
double timeFinal=6.653820100923719e+07;//Time period (s)
double deltaT;
//=1.330764020184744e+05;//Time step = time period/500 segments (s)
int numSteps = 50000;
deltaT = (timeFinal-timeInitial)/numSteps;
 
 // comparing completion time to MATLAB's time
 std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

// Our output function (yp)
  elements *yp = rk4sys(timeInitial,timeFinal,y0,deltaT);
  
  std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

  std::cout << "It took me " << time_span.count() << "seconds.";
  std::cout << std::endl;

  return 0;

// Output of r,theta,z,vr,vtheta,vz to a text file
  std::ofstream myfile;
  
  myfile.open ("orbitalMotion.txt");
  for(int i=0; i<numSteps;i++)
  {
  
  myfile << *(yp+i);
  }
  myfile.close();
  return 0;

delete(yp);
}