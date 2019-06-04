
#include <math.h>
#include "ode45.h"

elements calc_k(elements y, double h)
{
	
	elements k;
	k.r = calc_r(y);
	k.theta = calc_theta(y);
	k.z = calc_z(y);
	k.vr = calc_vr(y);
	k.vtheta = calc_vtheta(y);
	k.vz = calc_vz(y);
	return k;
}

elements calc_ymid(elements y, double h, elements k)
{
	elements ymid;
	/*
	ymid.r = y.r + k.r*h;
	ymid.theta = y.theta + k.theta*h;
	ymid.z = y.z + k.z*h;
	ymid.vr = y.vr + k.vr*h;
	ymid.vtheta = y.vtheta + k.theta*h;
	ymid.vz = y.vz + k.vz*h;
*/
	ymid = y+k*h;
	return ymid;
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