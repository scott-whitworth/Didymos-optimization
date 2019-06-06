#ifndef rk4sys_h
#define rk4sys_h
#include "ode45.h" // Utility functions for calc_k()


// Fourth-order Runge-Kutta algorthim for system of ODEs defined in ODE45.h
// Using linear time stepping 
    // If using non-linear time stepping, use two different functions
// Parameters: 
//      timeInitial: start time (s)
//      timeFinal: end time (s)
//      y0: initial conditions (position,velocity)
//      stepSize: time interval between data points (s)
// Output: A dynamic array of position and velocity sets, last entry is final conditions

// TODO: move allocation of output array outside of rk4sys()

template <class T> elements<T>* rk4sys(T timeInitial, T timeFinal, elements<T> y0,T stepSize, int numSteps, elements<T> *y);

#include "rk4sys.cpp"
#endif