#ifndef ode45_h
#define ode45_h

#include <ostream> // used in overload of stream output for elements
#include <iomanip> // setprecision(int)

#include "elements.h"
#include "coefficients.h"
#include "calcThrust.h"

#define AU 1.49597870691e11// units: m; used to convert meters to astronomical units
#define constG 6.67430e-11/(AU*AU*AU) //units: AU^3/(s^2 * kg); gravitational constant- used to calculate the gravitational force
#define massSun 1.98847e30//kg


//Calculates the corresponding k for the Runge-Kutta computation
// Units for k
//      k.r = au
//      k.theta = rad
//      k.z = au
//      k.vr = au/s
//      k.vtheta = rad/s
//      k.z = au/s
// Input:
//      y: current position and velocity conditions
//      h(time step): time interval between data points (s)
//      coeff: coefficients structure components for the current time stamp
//      accel: acceleration of the spacecraft (au/s^s)
//      curTime: current time stamp (s)
//      totalTime: the complete time frame of the simulation (s), used to normalize curTime
//Output: returns k1,k2,k3,k4 for y[n+1] calculation
template <class T> elements<T> calc_k(const T & h, const elements<T> & y, coefficients<T> & coeff, const T & accel, const T & curTime, const T & timeFinal);

// Dot = derivative of element with respect to time
// Utilities of calc_k(), calculates the element from current condition
// Parameter y: complete current condition

// Based on: y.vr
// Output: rDot
template <class T> T calcRate_r(const elements<T> & y);

// Based on: y.vtheta
// Output: thetaDot
template <class T> T calcRate_theta(const elements<T> & y);

// Based on: y.vz
// Output: zDot
template <class T> T calcRate_z(const elements<T> & y);

// other parameters for following equations:
//      coeff: coefficients structure components for the current time stamp
//      accel: acceleration of the spacecraft (au/s^s)
//      curTime: current time stamp (s)
//      totalTime: the complete time frame of the simulation (s), used to normalize curTime

// Based on: (-g * M_sun * r)  / (r^2 + z^2) ^ 3/2 + v_theta^2 / r + accel*cos(tau)*sin(gamma)
// Output: vrDot
template <class T> T calcRate_vr(const elements<T> & y, coefficients<T> & coeff, const T & accel, const T & curTime, const T & timeFinal);

// Based on: -vr*vtheta / r + accel*cos(tau)*cos(gamma)
// Output: vrDot
template <class T> T calcRate_vtheta(const elements<T> & y, coefficients<T> & coeff, const T & accel, const T & curTime, const T & timeFinal);

// Based on: (-g * M_sun * r)  / (r^2 + z^2) ^ 3/2 + + accel*sin(tau)
// Output: vrDot
template <class T> T calcRate_vz(const elements<T> & y, coefficients<T> & coeff, const T & accel, const T & curTime, const T & timeFinal);

#include "ode45.cpp"

#endif