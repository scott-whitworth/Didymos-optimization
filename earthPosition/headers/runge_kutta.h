#ifndef runge_kutta_h
#define runge_kutta_h
#include "motion_equations.h" // Utility functions for calc_k()


 

    // Comment on stepsize: Expected to be negative due to reverse integration
    // Output: writes in y the initial position of earth at the time of spacecraft launch based on an optimized trip time
    // Parameters: 
    //      timeInitial: start time (s)
    //      timeFinal: end time (s)
    //      times:an array that contains the time for each data point
    //      y0: initial conditions (position,velocity)
    //      stepSize: first time interval between data points (s)
    //      y: an array which contains the soultion to the dependent variable
    //      absTol: Sets the error tolerence for Runge-Kutta
    //      coeff: passes the structure containing the fourier coefficients for gamma and tau
    //      accel: spacecraft's acceleration (au/s^2)
    //      wetMass: mass of the spacecraft including fuel (kg)
    // To improve efficiency, the rk4 with single returns were split into two functions to avoid "if" statements, which are not prefered in CUDA.
template <class T> void rk4Reverse(const T & timeInitial, const T & timeFinal, const elements<T> & y0, 
T stepSize, elements<T> &y_new, const T & absTol);



// calculates k values 1 - 7 from equation and uses k values to find current and previous values of y
template <class T> void rkCalc(T & curTime, const T & timeFinal, T stepSize, elements<T> & y_new, elements<T> & error, elements<T> & k1,elements<T> & k2,
elements<T> & k3,elements<T> & k4, elements<T> & k5,elements<T> & k6,elements<T> & k7);


/**********************************************************************************************************************************/

// Calculates the scaling factor for the stepSize
// Parameters:
//      previous: The previous result of the Runge-Kutta 
//      difference: The new result minus the previous result (v-u)
//      absTol: Sets the error tolerence for Runge-Kutta
//      stepSize: time interval between data points (s)
// Output: Unitless scaling coefficient which changes the time step each iteration
template <class T> T calc_scalingFactor(const elements<T> & previous, const elements<T> & difference, const T & absTol, T & stepSize);

#include "runge_kutta.cpp"
#endif