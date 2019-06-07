#include "calcThrust.h"
#include <math.h>

template <class T> T calc_gamma(coefficients<T> const & coeff,T const & tt){
T gamma = coeff.gamma[0];

for(int i=1;i<=coeff.gamma.size()/2;i++){
gamma+=coeff.gamma[2*i-1]*cos(i*tt)+coeff.gamma[2*i]*sin(i*tt);
}
return gamma;
}

//----
template <class T> T calc_tau(coefficients<T> const & coeff,T const & tt){
T tau = coeff.tau[0];

for(int i=1;i<=coeff.tau.size()/2;i++){
tau+=coeff.tau[2*i-1]*cos(i*tt)+coeff.tau[2*i]*sin(i*tt);
}
return tau;
}

