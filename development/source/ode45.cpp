
#include <math.h>
#include "ode45.h"

elements calc_k(double h, elements y)
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

double calc_r(elements y)
{
	double r = y.vr;
	return r;
}

double calc_theta(elements y)
{
	double theta = y.vtheta / y.r;
	return theta;
}

double calc_z(elements y)
{
	double z = y.vz;
	return z;
}

double calc_vr(elements y)
{
	double vr=-constG * massSun * y.r / (pow(pow(y.r, 2) + pow(y.z, 2),(double)3/2)) + pow(y.vtheta,2) / y.r;
	return vr;
}

double calc_vtheta(elements y)
{
	double vtheta = -y.vr*y.vtheta / y.r;
	return vtheta;
}

double calc_vz(elements y)
{
	double vz = -constG * massSun * y.z / pow(pow(y.r, 2) + pow(y.z, 2),(double)3/2);
	return vz;
}