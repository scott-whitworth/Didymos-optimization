#ifndef rkParameters_h
#define rkParameters_h

#include "coefficients.h"
#include "elements.h"
#include "runge_kutta.h"

//struct to hold all the values required for the runge-kutta functions
template <class T> struct rkParameters {
    // Constructor which sets all the components
   __host__ __device__ rkParameters<T>(T timeFinal0, T accel0, T wetMass0, 
                   T r0, T theta0, T z0, T vr0, T vtheta0, T vz0,
                   T *gamma0, T *tau0, T *coast0, T coastThreshold0);

    __host__ __device__ rkParameters<T>();

    elements<T> y0;

    coefficients<T> coefficients;

    T timeFinal;
    T accel;
    T wetMass;

    //these are shared among all threads so they are passed seperately
    /*
    T timeInitial;
    T stepSize;
    T absTol;
    */

    // use the values in a rkParameters struct to call the serial rk4Simple() for comparison
    void parametersRK4Simple(T timeInitial, T stepSize, T absTol, elements<T> & y);
};

#include "rkParameters.cpp"
#endif