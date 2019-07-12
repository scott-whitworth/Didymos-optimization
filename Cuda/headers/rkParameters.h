#ifndef rkParameters_h
#define rkParameters_h

//structs
#include "coefficients.h"
#include "elements.h"

#include "runge_kutta.h"//used for rk4Simple()

//struct to hold all the values required for the runge-kutta functions
template <class T> struct rkParameters {
    // Constructor which sets all the components according to values taken in
   __host__ __device__ rkParameters<T>(T timeFinal0, T wetMass0, 
                   T r0, T theta0, T z0, T vr0, T vtheta0, T vz0, // elements<T>
                   T *gamma0, T *tau0, T *coast0, T coastThreshold0); // coefficients<T>

    // Alternate Constructor
    // timeFinal0           - Total time of mission (s)
    // wetMass              - Total mass of spacecraft at launch (kg)
    // initialCondition     - Initial position of the spacecraft at start of calculation (escape position/velocity)
    // coeff0               - Coefficients for the thrust angles and  acceleration
    __host__ __device__ rkParameters<T>(T timeFinal0, T wetMass0,
                  elements<T> initialCondition, coefficients<T> coeff0);          

    // constructor which sets everything to zero
    __host__ __device__ rkParameters<T>();

    // elements contains r, theta, z, vr, vtheta, and vz
    elements<T> y0;

    // coefficients contains arrays of Fourier series for:
    //    gamma
    //    tau
    //    coasting
    //and a value for coast_threshold
    coefficients<T> coeff;

    T timeFinal;
    T wetMass;

    // use the values in a rkParameters struct to call the serial rk4Simple() for comparison
    void parametersRK4Simple(T timeInitial, T stepSize, T absTol, elements<T> & y);
};

#include "rkParameters.cpp"
#endif