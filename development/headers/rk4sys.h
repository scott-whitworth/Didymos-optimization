#ifndef rk4sys_h
#define rk4sys_h
#include <vector>
#include "ode45.h"

//fourth-order Runge-Kutta for a system of ODEs
elements* rk4sys(double timeInitial, double timeFinal, elements y0, double stepSize);

void testKCalc(elements y0);

void testKCalc(elements y0);

#endif