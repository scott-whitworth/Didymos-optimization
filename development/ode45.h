#ifndef ode45_h
#define ode45_h


/* Finds the corresponding k for the Runge Kutta computation
 Input: the values of the vector y and the time step

 Output: return the k elements for the vector of equations
*/
elements k_calc(elements y, double h);

/* Calculates r, from the y vector

*/
double calc_r(elements y);


double calc_theta(elements y);

double calc_z(elements y);

double calc_vr(elements y);

double calc_vtheta(elements y);

double calc_vz(elements y);

#endif