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
template <class T> elements<T> calc_k(T const & h, elements<T> const & y, coefficients<T> const & coeff, T const & accel, T const & t, T const & timeFinal);

// Dot = derivative of element with respect to time
// Utilities of calc_k(), calculates the element from current condition
// Parameter y: complete current condition

// Based on: y.vr
// Output: rDot
template <class T> T calcRate_r(elements<T> const & y);

// Based on: y.vtheta
// Output: thetaDot
template <class T> T calcRate_theta(elements<T> const & y);

// Based on: y.vz
// Output: zDot
template <class T> T calcRate_z(elements<T> const & y);

// Based on: -constG * massSun * y.r / (pow(pow(y.r, 2) + pow(y.z, 2),(double)3/2)) + pow(y.vtheta,2) / y.r
// Output: vrDot
template <class T> T calcRate_vr(elements<T> const & y, coefficients<T> const & coeff, T const & accel, T const & t, T const & timeFinal);

// Based on: -y.vr*y.vtheta / y.r
// Output: vrDot
template <class T> T calcRate_vtheta(elements<T> const & y, coefficients<T> const & coeff, T const & accel, T const & t, T const & timeFinal);

// Based on: -constG * massSun * y.z / pow(pow(y.r, 2) + pow(y.z, 2),(double)3/2)
// Output: vrDot
template <class T> T calcRate_vz(elements<T> const & y, coefficients<T> const & coeff, T const & accel, T const & t, T const & timeFinal);

#include "ode45.cpp"

#endif