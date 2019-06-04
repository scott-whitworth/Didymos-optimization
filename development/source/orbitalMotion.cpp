#include "rk4sys.h"

int main()
{
elements y0;
y0.r = 3.150802646376772e+11;//Initial radial position of the asteroid (m)
y0.theta= -3.081519548404041;//Initial angular position of the asteroid (rad)
y0.z =  1.760293325286572e+10;//Initial off-plane position of the asteroid (m)
y0.vr = 4706.64912336045;//Initial radial velocity of the asteroid (m/s)
y0.vtheta= 16716.9055348804;//Intial azimuthal velocity of the asteroid (m/s)
y0.vz= -81.4453413932308;//Initial off-plane velocity of the asteroid (m/s)
double timeInitial=0;
double timeFinal=6.653820100923719e+07;
double h=1.330764020184744e+05;

elements *yp = rk4sys(timeInitial,timeFinal,y0,h);


// std::cout <<"hello!"<<std::endl;

// std::vector<elements> yp = rk4sys(tspan,y0,deltaT);

// std::cout << yp[0].r << yp[0].vr << yp[0].theta << yp.back().r;

testKCalc(y0);
}