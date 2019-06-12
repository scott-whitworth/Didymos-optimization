#ifndef coefficients_h
#define coefficients_h

// Utility to calculate the coefficients vector
//coefficients struct holds gamma and tau values

template <class T> struct coefficients {
    //setting the size of gamma
    //gamma will be evaluated as a fourier series
    const static int gammaSize=21;
    T gamma[gammaSize]; //in-plane coefficients angle
    
    //tau will be evaluated as a fourier series
    //setting the size of tau
    const static int tauSize=11;
    T tau[tauSize]; //out-of-plane coefficients angle

    //overload the stream output for elements used for writing to a file
    template <class U> friend std::ostream & operator<<(std::ostream & Str, const coefficients<T> & e); 
};

#include "coefficients.cpp"
#endif