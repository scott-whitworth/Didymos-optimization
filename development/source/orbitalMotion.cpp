#include "rk4sys.cpp"

int main()
{
elements y0;
y0.r = 3.150802646376772e+11;
y0.theta= -3.081519548404041;
y0.z =  1.760293325286572e+10;
y0.vr = 4706.64912336045;
y0.vtheta= 16716.9055348804;
y0.vz= -81.4453413932308;
std::vector <double> tspan;
tspan.front()=0;
tspan.back()=6.653820100923719e+07;
double h=1.330764020184744e+05;
//std::vector <double> *timeArray;

std::vector<elements> yp=rk4sys(tspan,y0,h);

}