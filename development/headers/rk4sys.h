#ifndef rk4sys_h
#define rk4sys_h
#include <vector>
#include "ode45.h"

//fourth-order RUnge-Kutta for a system of ODEs
std::vector<elements> rk4sys(std::vector<double> tspan, elements y0, double stepSize/*, std::vector<double> *tp*/);

#endif