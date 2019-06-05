#ifndef ode45_h
#define ode45_h

#include <ostream> // used in overload of stream output for elements
#include <iomanip> // setprecision(int)

#include "elements.h"

#define AU 1.49597870691e11// units: m; used to convert meters to astronomical units
#define constG 6.67430e-11/(AU*AU*AU) //units: AU^3/(s^2 * kg); gravitational constant- used to calculate the gravitational force
#define massSun 1.98847e30//kg

// TODO: Const ref changes

//Calculates the corresponding k for the Runge-Kutta computation
// Units for k
//      k.r = au
//      k.theta = rad
//      k.z = au
//      k.vr = au/s
//      k.vtheta = rad/s
//      k.z = au/s
// Input:
//       y: current position and velocity conditions
//       h(time step): time interval between data points (s)
//Output: returns k1,k2,k3,k4 for y[n+1] calculation
elements calc_k(double const & h, elements const & y);

// Dot = derivative of element with respect to time
// Utilities of calc_k(), calculates the element from current condition
// Parameter y: complete current condition

// Based on: y.vr
// Output: rDot
double calc_r(elements const & y);

// Based on: y.vtheta
// Output: thetaDot
double calc_theta(elements const & y);

// Based on: y.vz
// Output: zDot
double calc_z(elements const & y);

// Based on: -constG * massSun * y.r / (pow(pow(y.r, 2) + pow(y.z, 2),(double)3/2)) + pow(y.vtheta,2) / y.r
// Output: vrDot
double calc_vr(elements const & y);

// Based on: -y.vr*y.vtheta / y.r
// Output: vrDot
double calc_vtheta(elements const & y);

// Based on: -constG * massSun * y.z / pow(pow(y.r, 2) + pow(y.z, 2),(double)3/2)
// Output: vrDot
double calc_vz(elements const & y);

#endif