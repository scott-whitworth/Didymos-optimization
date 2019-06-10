#ifndef calcThrust_h
#define calcThrust_h

#include "coefficients.h"

template <class T> T calc_gamma(coefficients<T> const & coeff,T const & t, T const & timeFinal);

template <class T> T calc_tau(coefficients<T> const & coeff,T const & t, T const & timeFinal);

#include "calcThrust.cpp"
#endif