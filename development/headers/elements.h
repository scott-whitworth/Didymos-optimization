#ifndef elements_h
#define elements_h

#include <iostream>

//elements struct holds k values / dependent variable values in rk4sys
struct elements {
    //all in relation to the plane of the sun in cylindrical coordinates
    //Units are dependent upon context
    //positions
    double r; //radius (in plane)
    double theta; //angular position (in plane)
    double z; //axial position (out-of-plane)
    //velocities
    double vr; //radial velocity (in plane)
    double vtheta; //angular velocity (in plane)
    double vz; // axial velocity (out-of-plane)


    //overload operators to do math on all the elements in the struct seperately
    //Treating each element as a matrix operation

    //constructor which takes in an element
    elements operator+(const elements& e);
    elements operator-(const elements& e);
    elements operator*(const elements& e);
    elements operator/(const elements& e);

    //constructor which takes in an scalar
    elements operator*(const double& i);
    elements operator/(const double& i);

    //overload the stream output for elements used for writing to a file
    friend std::ostream & operator<<(std::ostream & Str, elements const & e); 
};

#endif