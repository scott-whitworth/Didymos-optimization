
#include <math.h>
#include "ode45.h"

elements calc_k(double h, elements y) //Calculation of k elements (k1,k2,k3,k4)
{
	
	elements k; 
	k.r = h*calc_r(y); // dimensional analysis (h=s, calc_r = au/s) = au
	k.theta = h*calc_theta(y); // dimensional analysis (h=s, calc_theta = rad/s) = rad
	k.z = h*calc_z(y);  // dimensional analysis (h=s, calc_z = au/s) = au
	k.vr = h*calc_vr(y); // dimensional analysis (h=s, calc_vr = au/s^2) = au/s
	k.vtheta = h*calc_vtheta(y); // dimensional analysis (h=s, calc_vtheta = rad/s^2) = rad/s
	k.vz = h*calc_vz(y); // dimensional analysis (h=s, calc_vr = au/s^2) = au/s
	return k;
}

/* Calculates rDot (dot = derivative of r with respect to time), from the y vector
 Input: the y vector
 Output : double rDot
*/
double calc_r(elements y)
{
	double r = y.vr;
	return r;
}

/* Calculates thetaDot, from the y vector
Input: the y vector
Output : double thetaDot
*/
double calc_theta(elements y)
{
	double theta = y.vtheta / y.r;
	return theta;
}

/* Calculates zDot, from the y vector
Input: the y vector
Output : double zDot
*/
double calc_z(elements y)
{
	double z = y.vz;
	return z;
}

/* Calculates vrDot, from the y vector
Input: the y vector
Output : double vrDot
*/
double calc_vr(elements y)
{
	double vr=-constG * massSun * y.r / (pow(pow(y.r, 2) + pow(y.z, 2),(double)3/2)) + pow(y.vtheta,2) / y.r;
	return vr;
}

/* Calculates vthetaDot, from the y vector
Input: the y vector
Output : double vthetaDot
*/
double calc_vtheta(elements y)
{
	double vtheta = -y.vr*y.vtheta / y.r;
	return vtheta;
}

/* Calculates vzDot, from the y vector
Input: the y vector
Output : double vzDot
*/
double calc_vz(elements y)
{
	double vz = -constG * massSun * y.z / pow(pow(y.r, 2) + pow(y.z, 2),(double)3/2);
	return vz;
}