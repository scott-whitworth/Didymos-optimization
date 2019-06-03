#include "ode45.h"

elements k_calc(elements y, double h)
{
	
	elements k[1];
	k.r = calc_r(y);
	k.theta = calc_theta(y);
	k.z = calc_z(y);
	return k1;
}

double calc_r(elements y)
{
	double r = y.vr;
	return r;
}

double calc_theta(elements y)
{
	theta = elements.vtheta / elements.r;
	return theta;
}

double calc_z(elements y)
{
	z = elements.vz;
	return z;
}

double calc_vr(elements y)
{
	vr=-constG * massSun * elements.r / (elements.r ^ 2 + elements.z ^ 2) ^ (3 / 2) + elements.vtheta ^ 2 / elements.r;
	return vr;
}

double calc_vtheta(elements y)
{
	vtheta = -elements.r*elements.vtheta / elements.r;
	return vtheta;
}

double calc_vz(elements y)
{
	vz = -constG * massSun * elements.z / (elements.r ^ 2 + elements.z ^ 2) ^ (3 / 2);
	return vz;
}
