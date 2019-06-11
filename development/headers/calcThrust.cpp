#include "calcThrust.h"
#include <math.h> // used for sine and cosine functions

template <class T> T calc_gamma(coefficients<T> const & coeff,T const & t, T const & timeFinal){
    T gamma = 0;
    gamma =coeff.gamma[0];

    for(int i=1;i<=(coeff.gammaSize-1)/2;i++){
        gamma+=coeff.gamma[2*i-1]*cos(i*t/timeFinal)+coeff.gamma[2*i]*sin(i*t/timeFinal);
    }
    return gamma;
}

//----
template <class T> T calc_tau(coefficients<T> const & coeff,T const & t, T const & timeFinal){
    T tau = 0;
    tau=coeff.tau[0];

    for(int i=1;i<=(coeff.tauSize-1)/2;i++){
        tau+=coeff.tau[2*i-1]*cos(i*t/timeFinal)+coeff.tau[2*i]*sin(i*t/timeFinal);
    }
    return tau;
}

//TODO: SC: These are the same function, they should be better functionalized.
//          If you are doing the same process, that is one function
//          I might suggest:
//          T calc_Series(const T & series[], const int & series_size, const T & curTime, const T & totalTime )
//          Then call this from calc_tau:
//          T calc_tau(const coefficients<T> & coeff, const T & t, const T & timeFinal){
//              return calc_Series(coeff.tau, coeff.tauSize,t,timeFinal);
//          }