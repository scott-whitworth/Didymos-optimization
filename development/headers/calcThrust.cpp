#include "calcThrust.h"
#include <math.h> // used for sine and cosine functions

template <class T> T calc_Series(T series[], const int series_size, const T & curTime, const T & totalTime){
    T coeff = series[0];

    for(int i=1;i<=(series_size-1)/2;i++){
        coeff+=series[2*i-1]*cos(i*curTime/totalTime)+series[2*i]*sin(i*curTime/totalTime);
    }
    return coeff;
}

template <class T> T calc_gamma(coefficients<T> & coeff,const T & t, const T & timeFinal){
   return calc_Series(coeff.gamma, coeff.gammaSize, t, timeFinal);
}

template <class T> T calc_tau(coefficients<T> & coeff, const T & t, const T & timeFinal){
    return calc_Series(coeff.tau, coeff.tauSize, t, timeFinal);
}
