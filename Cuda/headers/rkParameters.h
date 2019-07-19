#ifndef rkParameters_h
#define rkParameters_h

//structs
#include "coefficients.h"
#include "elements.h"

#include "runge_kutta.h"//used for rk4Simple()

//struct to hold all the values required for the runge-kutta functions
template <class T> struct rkParameters {
    
    //////////////////
    // Constructors //
    //////////////////
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
    
    /////////////
    // Members //
    /////////////
    // Initial Position/Velocity elements
    // Contains r, theta, z, vr, vtheta, and vz
    elements<T> y0;

    // Initial Optimization Coefficients
    // Contains arrays for Fourier series:
    //    gamma
    //    tau
    //    coasting
    // and a value for coast_threshold
    coefficients<T> coeff;

    // Final Time of simulation (s)
    T timeFinal;
    // Initial Wet Mass of spacecraft (kg)
    T wetMass;

    /////////////////////
    // Utility Methods //
    /////////////////////
    // Comparison function
    // Param: other - another rkParameter to be compared to
    // Param: comp_Thresh - comparison threshold
    // Returns true all elements of other are the same as *this, within the threshold comp_Thresh
    // Not exact, as there are a number of different magnitudes represented in parameters
    bool compare(const rkParameters<T> & other, T comp_Thresh);

    // Based on *this parameters, calculates the final position of the spacecraft using CPU based Methods
    // Param: timeInitial - start time of simulation (s)
    // Param: stepSize - first time interval between data points (s)
    // Param: absTol - error tolerence for Runge-Kutta
    // Output: y will contain the final position of the simulation, y initial value will be overwritten
    void parametersRK4Simple(T timeInitial, T stepSize, T absTol, elements<T> & y);
};

#include "rkParameters.cpp"
#endif