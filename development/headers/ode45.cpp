
#include <math.h>

#include "ode45.h"

template <class T> elements<T> calc_k(const T & h, const elements<T>  & y)
{
	elements<T> k; 
	k.r = h*calcRate_r(y); 
	k.theta = h*calcRate_theta(y); 
	k.z = h*calcRate_z(y);  
	k.vr = h*calcRate_vr(y); 
	k.vtheta = h*calcRate_vtheta(y); 
	k.vz = h*calcRate_vz(y); 
	return k;
}

template <class T> T calcRate_r(elements<T> const & y)
{
	return y.vr;
}

template <class T> T calcRate_theta(elements<T> const & y)
{
	return y.vtheta / y.r;
}

template <class T> T calcRate_z(elements<T> const & y)
{
	return y.vz;
}

template <class T> T calcRate_vr(elements<T> const & y)
{
	return -constG * massSun * y.r / (pow(pow(y.r, 2) + pow(y.z, 2),(T)3/2)) + pow(y.vtheta,2) / y.r;
}

template <class T> T calcRate_vtheta(elements<T> const & y)
{
	return -y.vr*y.vtheta / y.r;
}

template <class T> T calcRate_vz(elements<T> const & y)
{
	return -constG * massSun * y.z / pow(pow(y.r, 2) + pow(y.z, 2),(T)3/2);
}