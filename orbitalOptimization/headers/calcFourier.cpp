#include "calcFourier.h"
#include "acceleration.h"
#include <math.h> // used for sine and cosine functions
#include <iostream> // used for cout

template <class T> T calc_Series(T series[], const int series_size, const T & curTime, const T & timeFinal){
    T coeff = series[0];

    for(int i=1;i<=(series_size-1)/2;i++){
        coeff+=series[2*i-1]*cos(i*curTime/timeFinal)+series[2*i]*sin(i*curTime/timeFinal);
    }
    return coeff;
}

template <class T> T calc_gamma(coefficients<T> & coeff,const T & curTime, const T & timeFinal){
   return calc_Series(coeff.gamma, coeff.gammaSize, curTime, timeFinal);
}

template <class T> T calc_tau(coefficients<T> & coeff, const T & curTime, const T & timeFinal){
    return calc_Series(coeff.tau, coeff.tauSize, curTime, timeFinal);
}
