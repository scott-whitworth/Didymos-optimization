#ifndef ode45_h
#define ode45_h

#define constG 6.67430e-11
#define massSun 1.98847e30

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