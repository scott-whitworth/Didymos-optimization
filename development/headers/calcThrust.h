#ifndef calcThrust_h
#define calcThrust_h

#include "coefficients.h"

// Calculates Fourier series for a given input
// Parameters:
//         Series[]: specific array to be evaluated (ex: coeff.gamma)
//         series_size: size of the chosen array (ex: coeff.gammaSize)
//         curTime: current time (s) for calculated gamma
//         timeFinal: end time (s), used to normalize t
// output: Fourier series for a given input as a function of normalized time
template <class T> T calc_Series(T series[], const int series_size, const T & curTime, const T & timeFinal);

// Calculates gamma (in-plane angle) at a specific time using Fourier series function calc_Series
// Parameters:
//         coeff: coefficients structure, specifically the gamma components
//                coeff.gamma is an array of <T> values
//                coeff.gammaSize is the size of this array
//         t: current time (s) for calculated gamma
//         timeFinal: end time (s), used to normalize t
// output: in-plane angle derived from normalized time and gamma Fourier series
template <class T> T calc_gamma(const coefficients<T> & coeff,const T & curTime, const T & timeFinal);


// Calculates tau (out-of-plane angle) at a specific time using Fourier series function calc_Series
// Parameters:
//         coeff: coefficients structure, specifically the tau components
//                coeff.tau is an array of <T> values
//                coeff.tauSize is the size of this array
//         t: current time (s) for calculated tau
//         timeFinal: end time (s), used to normalize t
// output: in-plane angle derived from normalized time and tau Fourier series
template <class T> T calc_tau(const coefficients<T> & coeff,const T & curTime, const T & timeFinal);

#include "calcThrust.cpp"
#endif