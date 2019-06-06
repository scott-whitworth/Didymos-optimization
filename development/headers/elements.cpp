#include "elements.h"
#include <iomanip>


//overload operators to do math on all the elements in the struct seperately
//Treating each element as a matrix operation
//constructor which takes in an element
template <class T> 
elements<T> elements<T>::operator+(const elements<T>& e){
    elements<T> newElements;
    newElements.r = this->r + e.r;
    newElements.theta = this->theta + e.theta;
    newElements.z = this->z + e.z;
    newElements.vr = this->vr + e.vr;
    newElements.vtheta = this->vtheta + e.vtheta;
    newElements.vz = this->vz + e.vz;
    return newElements;
}

template <class T> 
elements<T> elements<T>::operator-(const elements& e){
    elements<T> newElements;
    newElements.r = this->r - e.r;
    newElements.theta = this->theta - e.theta;
    newElements.z = this->z - e.z;
    newElements.vr = this->vr - e.vr;
    newElements.vtheta = this->vtheta - e.vtheta;
    newElements.vz = this->vz - e.vz;
    return newElements;
}

template <class T> 
elements<T> elements<T>::operator*(const elements<T>& e){
    elements<T> newElements;
    newElements.r = this->r * e.r;
    newElements.theta = this->theta * e.theta;
    newElements.z = this->z * e.z;
    newElements.vr = this->vr * e.vr;
    newElements.vtheta = this->vtheta * e.vtheta;
    newElements.vz = this->vz * e.vz;
    return newElements;
}

template <class T> 
elements<T> elements<T>::operator/(const elements<T>& e){
    elements<T> newElements;
    newElements.r = this->r / e.r;
    newElements.theta = this->theta / e.theta;
    newElements.z = this->z / e.z;
    newElements.vr = this->vr / e.vr;
    newElements.vtheta = this->vtheta / e.vtheta;
    newElements.vz = this->vz / e.vz;
    return newElements;
}

//constructor which takes in an scalar
template <class T> 
elements<T> elements<T>::operator*(const T& i){
    elements<T> newElements;
    newElements.r = this->r * i;
    newElements.theta = this->theta * i;
    newElements.z = this->z * i;
    newElements.vr = this->vr * i;
    newElements.vtheta = this->vtheta * i;
    newElements.vz = this->vz * i;
    return newElements;
}

template <class T> elements<T> elements<T>::operator/(const T& i){
    elements<T> newElements;
    newElements.r = this->r / i;
    newElements.theta = this->theta / i;
    newElements.z = this->z / i;
    newElements.vr = this->vr / i;
    newElements.vtheta = this->vtheta / i;
    newElements.vz = this->vz / i;
    return newElements;
}

template <class T> std::ostream & operator<<(std::ostream & Str, elements<T> const & e) {
    Str << std::fixed;
    Str << std::setprecision(16); // number of decimals output into text file
    Str << e.r << "\t" << e.theta << "\t" << e.z << "\t" << e.vr << "\t" << e.vtheta << "\t" << e.vz << "\n";
    return Str;
}

