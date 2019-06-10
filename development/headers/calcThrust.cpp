#include "calcThrust.h"
#include <math.h>
#include <iostream>

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