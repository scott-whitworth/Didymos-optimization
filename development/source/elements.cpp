#include "elements.h"
#include <iomanip>


//overload operators to do math on all the elements in the struct seperately
//Treating each element as a matrix operation
//constructor which takes in an element
elements elements::operator+(const elements& e){
    elements newElements;
    newElements.r = this->r + e.r;
    newElements.theta = this->theta + e.theta;
    newElements.z = this->z + e.z;
    newElements.vr = this->vr + e.vr;
    newElements.vtheta = this->vtheta + e.vtheta;
    newElements.vz = this->vz + e.vz;
    return newElements;
}

    elements elements::operator-(const elements& e){
    elements newElements;
    newElements.r = this->r - e.r;
    newElements.theta = this->theta - e.theta;
    newElements.z = this->z - e.z;
    newElements.vr = this->vr - e.vr;
    newElements.vtheta = this->vtheta - e.vtheta;
    newElements.vz = this->vz - e.vz;
    return newElements;
}

elements elements::operator*(const elements& e){
    elements newElements;
    newElements.r = this->r * e.r;
    newElements.theta = this->theta * e.theta;
    newElements.z = this->z * e.z;
    newElements.vr = this->vr * e.vr;
    newElements.vtheta = this->vtheta * e.vtheta;
    newElements.vz = this->vz * e.vz;
    return newElements;
}

elements elements::operator/(const elements& e){
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
elements elements::operator*(const double& i){
    elements newElements;
    newElements.r = this->r * i;
    newElements.theta = this->theta * i;
    newElements.z = this->z * i;
    newElements.vr = this->vr * i;
    newElements.vtheta = this->vtheta * i;
    newElements.vz = this->vz * i;
    return newElements;
}

elements elements::operator/(const double& i){
    elements newElements;
    newElements.r = this->r / i;
    newElements.theta = this->theta / i;
    newElements.z = this->z / i;
    newElements.vr = this->vr / i;
    newElements.vtheta = this->vtheta / i;
    newElements.vz = this->vz / i;
    return newElements;
}

std::ostream & operator<<(std::ostream & Str, elements const & e) {
    Str << std::fixed;
    Str << std::setprecision(16); // number of decimals output into text file
    Str << e.r << "\t" << e.theta << "\t" << e.z << "\t" << e.vr << "\t" << e.vtheta << "\t" << e.vz << "\n";
    return Str;
}

