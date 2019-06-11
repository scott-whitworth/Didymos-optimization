#include "coefficients.h"
#include <iomanip> // setprecision(int)


//overload operators to do math on all the components in the struct seperately
//Treating each element as a matrix operation
//constructor which takes in an scalar
template <class T> 
coefficients<T> coefficients<T>::operator*(const T& i){
    coefficients<T> newCoefficients;

    //TODO: SC: Whoa! This is a huge issue. gamma is an array of T! You cannot treat gamma as a vector your are multiplying with a scalar like in matlab. This *will* cause errors
    //          A for loop is needed to multiply every value of gamma by i (same for tau)
    newCoefficients.gamma = this->gamma * i; 
    newCoefficients.tau = this->tau * i;
    return newCoefficients;
}

template <class T> coefficients<T> coefficients<T>::operator/(const T& i){
    coefficients<T> newCoefficients;
    newCoefficients.gamma = this->gamma / i;
    newCoefficients.tau = this->tau / i;
    return newCoefficients;
}

template <class T> std::ostream & operator<<(std::ostream & Str, coefficients<T> const & e) {
    Str << std::fixed;
    Str << std::setprecision(16); // number of decimals output into text file
    Str << e.gamma << "\t" << e.tau << "\n";
    return Str;
}