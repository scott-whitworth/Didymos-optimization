#ifndef elements_h
#define elements_h

#include <iostream> // used for cout

//elements struct holds k values / dependent variable values in rk4sys
template <class T> struct elements {
    // All in relation to the plane of the sun in cylindrical coordinates
    // Units are dependent upon context
    // Positions
    T r; //radius (in plane)
    T theta; //angular position (in plane)
    T z; //axial position (out-of-plane)
    //velocities
    T vr; //radial velocity (in plane)
    T vtheta; //angular velocity (in plane)
    T vz; // axial velocity (out-of-plane)

    // Constructor which sets all the components
    __host__ __device__ elements<T>(T r0, T theta0, T z0, T vr0, T vtheta0, T vz0);

    // Constructor which takes no values and sets everything to zero
    __host__ __device__ elements<T>();


    // Overload operators to do math on all the elements in the struct seperately
    // Treating each element as a matrix operation

    // Operator overloads, used in runge kutta for the math between elements
    __host__ __device__ elements<T> operator+(const elements<T>& e);
    __host__ __device__ elements<T> operator-(const elements<T>& e);
    __host__ __device__ elements<T> operator*(const elements<T>& e);
    __host__ __device__ elements<T> operator/(const elements<T>& e);

    // Operator overloads, for scalars
    __host__ __device__ elements<T> operator*(const T& i);
    __host__ __device__ elements<T> operator/(const T& i);

    // Overload the stream output for elements used for writing to a file
    template <class U> friend std::ostream & operator<<(std::ostream & Str, const elements<T> & e); 
};

#include "elements.cpp"

#endif