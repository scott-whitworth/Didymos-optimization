#ifndef elements_h
#define elements_h

#include <iostream>

//elements struct holds k values / dependent variable values in rk4sys
template <class T> struct elements {
    //all in relation to the plane of the sun in cylindrical coordinates
    //Units are dependent upon context
    //positions
    T r; //radius (in plane)
    T theta; //angular position (in plane)
    T z; //axial position (out-of-plane)
    //velocities
    T vr; //radial velocity (in plane)
    T vtheta; //angular velocity (in plane)
    T vz; // axial velocity (out-of-plane)


    //overload operators to do math on all the elements in the struct seperately
    //Treating each element as a matrix operation

    //constructor which takes in an element
    elements<T> operator+(const elements<T>& e);
    elements<T> operator-(const elements<T>& e);
    elements<T> operator*(const elements<T>& e);
    elements<T> operator/(const elements<T>& e);

    //constructor which takes in an scalar
    elements<T> operator*(const T& i);
    elements<T> operator/(const T& i);

    //overload the stream output for elements used for writing to a file
    template <class U> friend std::ostream & operator<<(std::ostream & Str, elements<T> const & e); 
};

#include "elements.cpp"

#endif