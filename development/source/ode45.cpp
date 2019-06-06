
#include <math.h>

#include "ode45.h"

elements calc_k(double const & h, elements const & y)
{
	elements k; 
	k.r = h*calcRate_r(y); 
	k.theta = h*calcRate_theta(y); 
	k.z = h*calcRate_z(y);  
	k.vr = h*calcRate_vr(y); 
	k.vtheta = h*calcRate_vtheta(y); 
	k.vz = h*calcRate_vz(y); 
	return k;
}

double calcRate_r(elements const & y)
{
	return y.vr;
}

double calcRate_theta(elements const & y)
{
	return y.vtheta / y.r;
}

double calcRate_z(elements const & y)
{
	return y.vz;
}

double calcRate_vr(elements const & y)
{
	return -constG * massSun * y.r / (pow(pow(y.r, 2) + pow(y.z, 2),(double)3/2)) + pow(y.vtheta,2) / y.r;
}

double calcRate_vtheta(elements const & y)
{
	return -y.vr*y.vtheta / y.r;
}

double calcRate_vz(elements const & y)
{
	return -constG * massSun * y.z / pow(pow(y.r, 2) + pow(y.z, 2),(double)3/2);
}