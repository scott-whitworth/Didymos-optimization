
#include <math.h>

#include "ode45.h"

elements calc_k(double const & h, elements const & y)
{
	elements k; 
	k.r = h*calc_r(y); 
	k.theta = h*calc_theta(y); 
	k.z = h*calc_z(y);  
	k.vr = h*calc_vr(y); 
	k.vtheta = h*calc_vtheta(y); 
	k.vz = h*calc_vz(y); 
	return k;
}

double calc_r(elements const & y)
{
	return y.vr;
}

double calc_theta(elements const & y)
{
	return y.vtheta / y.r;
}

double calc_z(elements const & y)
{
	return y.vz;
}

double calc_vr(elements const & y)
{
	return -constG * massSun * y.r / (pow(pow(y.r, 2) + pow(y.z, 2),(double)3/2)) + pow(y.vtheta,2) / y.r;
}

double calc_vtheta(elements const & y)
{
	return -y.vr*y.vtheta / y.r;
}

double calc_vz(elements const & y)
{
	return -constG * massSun * y.z / pow(pow(y.r, 2) + pow(y.z, 2),(double)3/2);
}