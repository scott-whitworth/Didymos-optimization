#ifndef coefficients_h
#define coefficients_h

// Utility to calculate the coefficients vector
//coefficients struct holds acceleration, gamma, and tau values
//TODO: SC: This could use a lot of clarifying comments. At no point do you reference these being a fourier set
template <class T> struct coefficients {
    //setting the size of gamma
    const static int gammaSize=21;
    T gamma[gammaSize]; //in-plane coefficients angle
    //setting the size of tau
    const static int tauSize=11;
    T tau[tauSize]; //out-of-plane coefficients angle

    //overload the stream output for elements used for writing to a file
    template <class U> friend std::ostream & operator<<(std::ostream & Str, const coefficients<T> & e); 
};

#include "coefficients.cpp"
#endif