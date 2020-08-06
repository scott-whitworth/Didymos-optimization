#ifndef rk4sys_h
#define rk4sys_h
#include "ode45.h" // Utility functions for calc_k()

//TODO: SC: This needs to be updated to include information on coeff, gamma and tau

// Fourth-order Runge-Kutta algorthim for system of ODEs defined in ODE45.h
// Using adaptive time stepping  
// Parameters: 
//      timeInitial: start time (s)
//      timeFinal: end time (s)
//      times:an array that contains the time for each data point
//      y0: initial conditions (position,velocity)
//      stepSize: time interval between data points (s)
//      y: an array which contains the soultion to the dependent variable
//      absTol: Sets the error tolerence for Runge-Kutta
// Output: A dynamic array of position and velocity sets, last entry is final conditions
template <class T> void rk4sys(const T & timeInitial, const T & timeFinal, T *times, const elements<T> & y0,  T stepSize, elements<T> *y, const T & absTol, coefficients<T> coeff, const T & accel, T *gamma, T *tau);

//take a time step in the runge-kutta algorithm
//called repeatedly in rk4sys until the final time is reached
template <class T> void rkStep(T *curTime, const T & timeInitial, const T & timeFinal, T *times, T stepSize, elements<T> *y, const T & absTol, coefficients<T> & coeff, const T & accel, T *gamma, T *tau, int *iteration);

//calculate the runge-kutta equations for the current time
//called from rkStep
template <class T> void rkCalc(T *curTime, const T & timeFinal, T stepSize, elements<T> *y, coefficients<T> & coeff, const T & accel, elements<T> & v, elements<T> & u, int iteration);

// Calculates the scaling factor for the stepSize
// Parameters:
//      previous: The previous result of the Runge-Kutta 
//      difference: The new result minus the previous result (v-u)
//      absTol: Sets the error tolerence for Runge-Kutta
//      stepSize: time interval between data points (s)
// Output: Unitless scaling coefficient which changes the time step each iteration
template <class T> T calc_scalingFactor(const elements<T> & previous , const elements<T> & difference, const T & absTol, T & stepSize);

#include "rk4sys.cpp"
#endif