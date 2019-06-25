#ifndef coefficients_h
#define coefficients_h

//Utility to calculate the coefficients vector
//Coefficients struct holds gamma and tau values
template <class T> struct coefficients {
    
    // fourth order fourier series
    // setting the size of gamma (in-plane coefficients angle)
    const static int gammaSize=9; 
    T gamma[gammaSize]; 

    // second order fourier series
    // setting the size of tau (out-of-plane coefficients angle)
    const static int tauSize=3; 
    T tau[tauSize]; 

    const static int coastSize=5;
    T coast[coastSize];

    T coastThreshold;

    //overload the stream output for elements used for writing to a file
    template <class U> friend std::ostream & operator<<(std::ostream & Str, const coefficients<T> & e); 
};

#include "coefficients.cpp"
#endif