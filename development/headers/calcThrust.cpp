#include "calcThrust.h"
#include <math.h>
#include <iostream>

// ask scott how to do this with out hard coding in 10 and 5

template <class T> T calc_gamma(coefficients<T> const & coeff,T const & t){
T gamma = 0;
gamma =coeff.gamma[0];

for(int i=1;i<=10;i++){
gamma+=coeff.gamma[2*i-1]*cos(i*t)+coeff.gamma[2*i]*sin(i*t);
}
std::cout<<gamma<<'\n';
return gamma;
}

//----
template <class T> T calc_tau(coefficients<T> const & coeff,T const & t){
T tau = 0;
tau=coeff.tau[0];

for(int i=1;i<=5;i++){
tau+=coeff.tau[2*i-1]*cos(i*t)+coeff.tau[2*i]*sin(i*t);
}
std::cout<<tau<<'\n';
return tau;
}