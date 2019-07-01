
#include <math.h> // used for sine, cosine, and pow functions
#include "motion_equations.h"

template <class T> elements<T> calc_k(const T & h, const elements<T>  & y, coefficients<T> & coeff, const T & accel, const T & curTime, const T & timeFinal)
{
	return elements<T>( h*calcRate_r(y), h*calcRate_theta(y), h*calcRate_z(y), 
	h*calcRate_vr(y,coeff,accel,curTime, timeFinal), h*calcRate_vtheta(y,coeff,accel,curTime, timeFinal),  h*calcRate_vz(y,coeff,accel,curTime, timeFinal));
}

template <class T> T calcRate_r(const elements<T> & y)
{
	return y.vr;
}

template <class T> T calcRate_theta(const elements<T> & y)
{
	return y.vtheta / y.r;
}

template <class T> T calcRate_z(const elements<T> & y)
{
	return y.vz;
}

template <class T> T calcRate_vr(const elements<T> & y, coefficients<T> & coeff, const T & accel, const T & curTime, const T & timeFinal)
{
	return (-constG * massSun * y.r / (pow(pow(y.r, 2) + pow(y.z, 2), (T)3/2))) + (pow(y.vtheta,2) / y.r) +
	 (accel*cos(calc_tau(coeff,curTime, timeFinal))*sin(calc_gamma(coeff,curTime, timeFinal)));
	
}

template <class T> T calcRate_vtheta(const elements<T> & y, coefficients<T> & coeff, const T & accel, const T & curTime, const T & timeFinal)
{
	return -y.vr*y.vtheta / y.r + accel*cos(calc_tau(coeff,curTime, timeFinal))*cos(calc_gamma(coeff,curTime, timeFinal));

}

template <class T> T calcRate_vz(const elements<T> & y, coefficients<T> & coeff, const T & accel, const T & curTime, const T & timeFinal)
{
	return (-constG * massSun * y.z / pow(pow(y.r, 2) + pow(y.z, 2), (T)3/2)) + accel*sin(calc_tau(coeff,curTime, timeFinal));
	
}