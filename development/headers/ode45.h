#ifndef ode45_h
#define ode45_h

#include <ostream> // used in overload of stream output for elements
#include <iomanip> // setprecision(int)

#define AU 1.49597870691e11// units: m; used to convert meters to astronomical units
#define constG 6.67430e-11/(AU*AU*AU) //units: AU^3/(s^2 * kg); gravitational constant- used to calculate the gravitational force
#define massSun 1.98847e30//kg

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
    elements operator+(const elements& e){
        elements newElements;
        newElements.r = this->r + e.r;
        newElements.theta = this->theta + e.theta;
        newElements.z = this->z + e.z;
        newElements.vr = this->vr + e.vr;
        newElements.vtheta = this->vtheta + e.vtheta;
        newElements.vz = this->vz + e.vz;
        return newElements;
    }

     elements operator-(const elements& e){
        elements newElements;
        newElements.r = this->r - e.r;
        newElements.theta = this->theta - e.theta;
        newElements.z = this->z - e.z;
        newElements.vr = this->vr - e.vr;
        newElements.vtheta = this->vtheta - e.vtheta;
        newElements.vz = this->vz - e.vz;
        return newElements;
    }

    elements operator*(const elements& e){
        elements newElements;
        newElements.r = this->r * e.r;
        newElements.theta = this->theta * e.theta;
        newElements.z = this->z * e.z;
        newElements.vr = this->vr * e.vr;
        newElements.vtheta = this->vtheta * e.vtheta;
        newElements.vz = this->vz * e.vz;
        return newElements;
    }

    elements operator/(const elements& e){
        elements newElements;
        newElements.r = this->r / e.r;
        newElements.theta = this->theta / e.theta;
        newElements.z = this->z / e.z;
        newElements.vr = this->vr / e.vr;
        newElements.vtheta = this->vtheta / e.vtheta;
        newElements.vz = this->vz / e.vz;
        return newElements;
    }

    //constructor which takes in an scalar
    elements operator*(const double& i){
        elements newElements;
        newElements.r = this->r * i;
        newElements.theta = this->theta * i;
        newElements.z = this->z * i;
        newElements.vr = this->vr * i;
        newElements.vtheta = this->vtheta * i;
        newElements.vz = this->vz * i;
        return newElements;
    }

    elements operator/(const double& i){
        elements newElements;
        newElements.r = this->r / i;
        newElements.theta = this->theta / i;
        newElements.z = this->z / i;
        newElements.vr = this->vr / i;
        newElements.vtheta = this->vtheta / i;
        newElements.vz = this->vz / i;
        return newElements;
    }
};

//overload the stream output for elements used for writing to a file
inline std::ostream & operator<<(std::ostream & Str, elements const & e) {
    Str << std::fixed;
    Str << std::setprecision(16); // number of decimals output into text file
    Str << e.r << "\t" << e.theta << "\t" << e.z << "\t" << e.vr << "\t" << e.vtheta << "\t" << e.vz << "\n";
    return Str;
}

// TODO: Pull apart elements definition from ODE45 functions into a new header and cpp
// TODO: Const ref changes

// Units for k
// k.r = au
// k.theta = rad
// k.z = au
// k.vr = au/s
// k.vtheta = rad/s
// k.z = au/s

/* Calculates the corresponding k for the Runge Kutta computation
 Input: the values of the vector y (r,theta,z,vr,vtheta,vz) and the time step
 Output: return the k elements for the vector of equations (k1,k2,k3,k4)
*/
elements calc_k(double h, elements y);

// Dot = derivative of element with respect to time

// TODO: redo all other comments based on r_calc() and modify k_calc() comments

// Utility of calc_k(), calculates the rDot element from current condition
// Based on: vr
// Parameter y: complete current condition
// Output: rDot
double calc_r(elements y);

/* Calculates thetaDot, from the y vector
Input: the y vector
Output : double thetaDot
*/
double calc_theta(elements y);

/* Calculates zDot, from the y vector
Input: the y vector
Output : double zDot
*/
double calc_z(elements y);

/* Calculates vrDot, from the y vector
Input: the y vector
Output : double vrDot
*/
double calc_vr(elements y);

/* Calculates vthetaDot, from the y vector
Input: the y vector
Output : double vthetaDot
*/
double calc_vtheta(elements y);

/* Calculates vzDot, from the y vector
Input: the y vector
Output : double vzDot
*/
double calc_vz(elements y);

#endif