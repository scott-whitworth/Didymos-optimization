#include "rkParameters.h"

template <class T> rkParameters<T>::rkParameters(T timeFinal0, T accel0, T wetMass0, 
                                                 T r0, T theta0, T z0, T vr0, T vtheta0, T vz0, // elements<T>
                                                 T gamma0[], T tau0[], T coast0[], T coastThreshold0) // coefficients<T>
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

    coefficients.gamma = gamma0;
    coefficients.tau = tau0;
    coefficients.coast = coast0;
    coefficients.coastThreshold = coastThreshold0;
}