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
	return k1;
}

elements calc_ymid(elements y, double h, elements k)
{
	elements ymid;
	ymid.r = y.r + k.r*h;
	ymid.theta = y.theta + k.theta*h;
	ymid.z = y.z + k.z*h;
	ymid.vr = y.vr + k.vr*h;
	ymid.vtheta = y.vtheta + k.theta*h;
	ymid.vz = y.vz + k.vz*h;

	return ymid;
}

double calc_r(elements y)
{
	r = elements.vr;
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
