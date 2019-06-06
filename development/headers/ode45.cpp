
#include <math.h>

#include "ode45.h"

template <class T> elements<T> calc_k(const T & h, const elements<T>  & y,T accel,T tau, T gamma)
{
	elements<T> k; 
	k.r = h*calcRate_r(y); 
	k.theta = h*calcRate_theta(y); 
	k.z = h*calcRate_z(y);  
	k.vr = h*calcRate_vr(y,accel,tau,gamma); 
	k.vtheta = h*calcRate_vtheta(y,accel,tau,gamma); 
	k.vz = h*calcRate_vz(y,accel,tau); 
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

template <class T> T calcRate_vr(elements<T> const & y, T accel, T tau, T gamma)
{
	return (-constG * massSun * y.r / (pow(pow(y.r, 2) + pow(y.z, 2), (T)3/2))) + (pow(y.vtheta,2) / y.r) + (accel*cos(tau)*sin(gamma));
}

template <class T> T calcRate_vtheta(elements<T> const & y, T accel, T tau, T gamma)
{
	return -y.vr*y.vtheta / y.r + accel*cos(tau)*cos(gamma);
}

template <class T> T calcRate_vz(elements<T> const & y, T accel, T tau)
{
	return (-constG * massSun * y.z / pow(pow(y.r, 2) + pow(y.z, 2), (T)3/2)) + accel*sin(tau);
}