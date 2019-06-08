#ifndef coefficients_h
#define coefficients_h

// Utility to calculate the coefficients vector
//coefficients struct holds acceleration, gamma, and tau values
template <class T> struct coefficients {
    T gamma[21]; //in-plane coefficients angle
    T tau[11]; //out-of-plane coefficients angle

    //overload operators to do math on all the components in the struct seperately
    //constructor which takes in an scalar
    //multiplies each component of coefficients by a scalar value
    coefficients<T> operator*(const T& i);
    //divides each component of coefficients by a scalar value
    coefficients<T> operator/(const T& i);

    //overload the stream output for elements used for writing to a file
    template <class U> friend std::ostream & operator<<(std::ostream & Str, coefficients<T> const & e); 
};

#include "coefficients.cpp"
#endif