#ifndef rkParameters_h
#define rkParameters_h

#include "coefficients.h"
#include "elements.h"

//struct to hold all the values required for the runge-kutta functions
template <class T> struct rkParameters {
    // Constructor which sets all the components
   rkParameters<T>(T timeFinal0, T stepSize0, T absTol0, T accel0, T wetMass0, 
                   T r0, T theta0, T z0, T vr0, T vtheta0, T vz0,
                   T gamma0[], T tau0[], T coast0[], T coastThreshold0);



    elements<T> y0, y;

    coefficients<T> coefficients;

    T timeFinal;

    T stepSize;

    T absTol;

    T accel;

    T wetMass;
};

#endif