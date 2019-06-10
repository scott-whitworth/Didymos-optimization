#ifndef rk4sys_h
#define rk4sys_h
#include "ode45.h" // Utility functions for calc_k()

// TODO: move allocation of output array outside of rk4sys()

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
template <class T> elements<T>* rk4sys(T timeInitial, T timeFinal,T *times, elements<T> y0, T stepSize, elements<T> *y, T absTol, coefficients<T> coeff, T accel,T *gamma, T *tau);

// Calculates the scaling factor for the stepSize
// Parameters:
//      previous: The previous result of the Runge-Kutta 
//      difference: The new result minus the previous result (v-u)
//      absTol: Sets the error tolerence for Runge-Kutta
//      stepSize: time interval between data points (s)
// Output: Unitless scaling coefficient which changes the time step each iteration
template <class T> T calc_scalingFactor(elements<T> previous ,elements<T> difference, T absTol, T stepSize);

#include "rk4sys.cpp"
#endif