// Didymos-Optimization_Project:
// Last Editor: Lauren and Ben
// Tasks Completed: 
    // Created calc_coast()

#include "calcFourier.h"
#include <math.h> // Used for sine and cosine functions
#include <iostream> // Used for cout

template <class T> T calc_Series(const T series[], const int series_size, const T & curTime, const T & timeFinal){
    T coeff = series[0];
    T curTimeRatio = curTime / timeFinal;

    // f(x) = a_0 + sum{a_n*cos(2*pi*n*t/L)+b_n*sin(2*pi*n*t/L)}
    for(int i=1;i<=(series_size-1)/2;i++)
    {
        coeff+=series[2*i-1]*cos(2*M_PI*i*curTimeRatio)+series[2*i]*sin(2*M_PI*i*curTimeRatio);
    }
    return coeff;
}

template <class T> T calc_gamma(coefficients<T> & coeff,const T & curTime, const T & timeFinal){
   return calc_Series(coeff.gamma, coeff.gammaSize, curTime, timeFinal);
}

template <class T> T calc_tau(coefficients<T> & coeff, const T & curTime, const T & timeFinal){
    return calc_Series(coeff.tau, coeff.tauSize, curTime, timeFinal);
}

template <class T> bool calc_coast(coefficients<T> & coeff, const T & curTime, const T & timeFinal)
{
    // Use the fourier series for the coasting coefficients, then take the cos^2(coasting)
    T coastValue = pow(sin(calc_Series(coeff.coast, coeff.coastSize, curTime, timeFinal)),2);
    // If coastValue is above the optimized threshold we return true for not coasting
    if(coastValue>=coeff.coastThreshold){
        return true;
    }
    // Otherwise false
    else{
        return false;
    }
}