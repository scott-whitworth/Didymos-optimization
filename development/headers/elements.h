#ifndef elements_h
#define elements_h

// TODO: Pull apart elements definition into cpp

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

#endif