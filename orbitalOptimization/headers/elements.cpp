// Didymos-Optimization_Project:
// Last Editor: Lauren
// Tasks Completed: 
    // Changed constructors to take in const values

#include "elements.h"
#include <iomanip> // setprecision(int)


// Constrctors

// Sets starting values as given
template <class T>
elements<T>::elements(T r0, T theta0, T z0, T vr0, T vtheta0, T vz0){
    r = r0;
    theta = theta0;
    z = z0;
    vr = vr0;
    vtheta = vtheta0;
    vz = vz0;
}

// Sets values to default of zero
template <class T>
elements<T>::elements(){
    r = static_cast<T>(0);
    theta = static_cast<T>(0);
    z = static_cast<T>(0);
    vr = static_cast<T>(0);
    vtheta = static_cast<T>(0);
    vz = static_cast<T>(0);
}


// Overload operators to do math on all the elements in the struct seperately WITH CONST
// Treating each element as a matrix operation

// Constructor which takes in an element
template <class T> 
elements<T> elements<T>::operator+(const elements<T> & e)const{
    return elements<T>(this->r + e.r, this->theta + e.theta, this->z + e.z, this->vr + e.vr, this->vtheta + e.vtheta, this->vz + e.vz);
}

template <class T> 
elements<T> elements<T>::operator-(const elements & e)const{
    return elements<T>(this->r - e.r, this->theta - e.theta, this->z - e.z, this->vr - e.vr, this->vtheta - e.vtheta, this->vz - e.vz);
}

template <class T> 
elements<T> elements<T>::operator*(const elements<T> & e)const{
    return elements<T>(this->r * e.r, this->theta * e.theta, this->z * e.z, this->vr * e.vr, this->vtheta * e.vtheta, this->vz * e.vz);
}

template <class T> 
elements<T> elements<T>::operator/(const elements<T> & e)const{
    return elements<T>(this->r / e.r, this->theta / e.theta, this->z / e.z, this->vr / e.vr, this->vtheta / e.vtheta, this->vz / e.vz);
}

// Constructor which takes in an scalar
template <class T> 
elements<T> elements<T>::operator*(const T & i)const{
    return elements<T>(this->r * i,  this->theta * i, this->z * i, this->vr * i, this->vtheta * i, this->vz * i);
}

template <class T> elements<T> elements<T>::operator/(const T & i)const{
    return elements<T>( this->r / i, this->theta / i, this->z / i, this->vr / i, this->vtheta / i, this->vz / i);
}

template <class T> std::ostream & operator<<(std::ostream & Str, const elements<T> & e) {
    Str << std::fixed;
    Str << std::setprecision(15); // Number of decimals output into text file
    Str << std::scientific;
    Str << e.r << "\t" << e.theta << "\t" << e.z << "\t" << e.vr << "\t" << e.vtheta << "\t" << e.vz << "\n";
    return Str;
}

