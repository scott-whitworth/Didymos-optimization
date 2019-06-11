
#include <math.h> // used for sine, cosine, and pow functions
#include "ode45.h"

template <class T> elements<T> calc_k(const T & h, const elements<T>  & y, const coefficients<T> & coeff, const T & accel,  T const  & t, T const & timeFinal)
{
	return elements<T>( h*calcRate_r(y), h*calcRate_theta(y), h*calcRate_z(y), 
	h*calcRate_vr(y,coeff,accel,t, timeFinal), h*calcRate_vtheta(y,coeff,accel,t, timeFinal),  h*calcRate_vz(y,coeff,accel,t, timeFinal));
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

template <class T> T calcRate_vr(elements<T> const & y,  coefficients<T> const & coeff, T const & accel,T const & t, T const & timeFinal)
{
	return (-constG * massSun * y.r / (pow(pow(y.r, 2) + pow(y.z, 2), (T)3/2))) + (pow(y.vtheta,2) / y.r) + (accel*cos(calc_tau(coeff,t, timeFinal))*sin(calc_gamma(coeff,t, timeFinal)));
}

template <class T> T calcRate_vtheta(elements<T> const & y, coefficients<T> const & coeff, T const & accel, T const & t, T const & timeFinal)
{
	return -y.vr*y.vtheta / y.r + accel*cos(calc_tau(coeff,t, timeFinal))*cos(calc_gamma(coeff,t, timeFinal));
}

template <class T> T calcRate_vz(elements<T> const & y,coefficients<T> const & coeff, T const & accel, T const & t, T const & timeFinal)
{
	return (-constG * massSun * y.z / pow(pow(y.r, 2) + pow(y.z, 2), (T)3/2)) + accel*sin(calc_tau(coeff,t, timeFinal));
}