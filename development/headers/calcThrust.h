#ifndef calcThrust_h
#define calcThrust_h

#include "coefficients.h"
#include <cmath>


template <class T> T calc_gamma(coefficients<T> const & coeff,T const & tt);

template <class T> T calc_tau(coefficients<T> const & coeff,T const & tt);

#endif