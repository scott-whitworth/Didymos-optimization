#ifndef ode45_h
#define ode45_h
#define constG 6.67430e-11
#define massSun 1.98847e30


//elements struct holds k values / dependent variable values in rk4sys
struct elements {
    double r;
    double theta;
    double z;
    double vr;
    double vtheta;
    double vz;

    //overload operators to do math on all the elements in the struct seperately
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

    elements operator*(const double i){
        elements newElements;
        newElements.r = this->r * i;
        newElements.theta = this->theta * i;
        newElements.z = this->z * i;
        newElements.vr = this->vr * i;
        newElements.vtheta = this->vtheta * i;
        newElements.vz = this->vz * i;
        return newElements;
    }

    elements operator/(const double i){
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



/* Finds the corresponding k for the Runge Kutta computation
 Input: the values of the vector y and the time step

 Output: return the k elements for the vector of equations
*/
elements calc_k(elements y, double h);

/*
*/
elements calc_ymid(elements y, double h, elements k);

/* Calculates r, from the y vector
 Input: the y vector
 Output : double r
*/
double calc_r(elements y);

/* Calculates theta, from the y vector
Input: the y vector
Output : double theta
*/
double calc_theta(elements y);

/* Calculates z, from the y vector
Input: the y vector
Output : double z
*/
double calc_z(elements y);

/* Calculates vr, from the y vector
Input: the y vector
Output : double vr
*/
double calc_vr(elements y);

/* Calculates vtheta, from the y vector
Input: the y vector
Output : double vtheta
*/
double calc_vtheta(elements y);

/* Calculates vz, from the y vector
Input: the y vector
Output : double vz
*/
double calc_vz(elements y);

#endif