#include "rkParameters.h"

template <class T> rkParameters<T>::rkParameters(T timeFinal0, T accel0, T wetMass0, 
                                                 T r0, T theta0, T z0, T vr0, T vtheta0, T vz0, // elements<T>
                                                 T *gamma0, T *tau0, T *coast0, T coastThreshold0) // coefficients<T>
{
    timeFinal = timeFinal0;
    accel = accel0;
    wetMass = wetMass0;

    y0.r = r0;
    y0.theta = theta0;
    y0.z = z0;
    y0.vr = vr0;
    y0.vtheta = vtheta0;
    y0.vz = vz0;

    for(int i = 0; i < coefficients.gammaSize; i++){
        coefficients.gamma[i] = gamma0[i];
    }
    for(int i = 0; i < coefficients.tauSize; i++){
        coefficients.tau[i] = tau0[i];
    }
    for(int i = 0; i < coefficients.coastSize; i++){
        coefficients.coast[i] = coast0[i];
    }

    coefficients.coastThreshold = coastThreshold0;
}

template <class T> rkParameters<T>::rkParameters()
{
    timeFinal = 0;
    accel = 0;
    wetMass = 0;

    y0.r = 0;
    y0.theta = 0;
    y0.z = 0;
    y0.vr = 0;
    y0.vtheta = 0;
    y0.vz = 0;

    coefficients.gamma = NULL;
    coefficients.tau = NULL;
    coefficients.coast = NULL;
    coefficients.coastThreshold = 0;
}

template <class T> void rkParameters<T>::parametersRK4Simple(T timeInitial, T stepSize, T absTol, elements<T> & y){
    rk4Simple(timeInitial, timeFinal, y0, stepSize, y, absTol, coefficients, accel, wetMass);
}