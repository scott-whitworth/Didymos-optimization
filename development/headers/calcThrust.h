#ifndef calcThrust_h
#define calcThrust_h

#include "coefficients.h"


// Fourier series used to calculate gamma (in-plane angle)
// input:
//         coeff.gamma: coefficients structure component, used to set initial gamma
//         t: time stamp for a given gamma (s)
//         timeFinal: end time (s), used to normalize t
// output: in-plane angle as a function of normalized time
template <class T> T calc_gamma(coefficients<T> const & coeff,T const & t, T const & timeFinal);


// Fourier series used to calculate tau (out-of-plane angle)
// input:
//         coeff.tau: coefficients structure component, used to set initial tau
//         t: time stamp for a given tau (s)
//         timeFinal: end time (s), used to normalize t
// output: out-of-plane angle as a function of normalized time
template <class T> T calc_tau(coefficients<T> const & coeff,T const & t, T const & timeFinal);

#include "calcThrust.cpp"
#endif