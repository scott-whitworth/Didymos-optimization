// Didymos-Optimization_Project:
// Last Editor: Mateo and Lauren
// Tasks Completed: 
	// No recent changes
	
#include <math.h> // Used for sine, cosine, and pow functions
#include "motion_equations.h"

//////////////////////////////////////////////////////////////
// "Free motion" equations due to thrust
//////////////////////////////////////////////////////////////

template <class T> elements<T> calc_k(const T & h, const elements<T>  & y, coefficients<T> & coeff, const T & accel, const T & curTime, const T & timeFinal)
{
	return elements<T>( h*calcRate_r(y), h*calcRate_theta(y), h*calcRate_z(y), 
	h*calcRate_vr(y,coeff,accel,curTime, timeFinal), h*calcRate_vtheta(y,coeff,accel,curTime, timeFinal),  h*calcRate_vz(y,coeff,accel,curTime, timeFinal));
}

// Calc k for cases with no thrust
// Changed order to match the order in GPU code
template <class T> elements<T> calc_k_earth(const T & h, const elements<T>  & y, const T & curTime, const T & timeFinal)
{
	return elements<T>( h*calcRate_r(y), h*calcRate_theta(y), h*calcRate_z(y), 
	h*calcRate_vr_earth(y,curTime, timeFinal), h*calcRate_vtheta_earth(y,curTime, timeFinal),  h*calcRate_vz_earth(y,curTime, timeFinal));
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
	return (-constG * massSun * y.r / (pow(pow(y.r, 2) + pow(y.z, 2), 1.5))) + (pow(y.vtheta,2) / y.r) +
	 (accel*cos(calc_tau(coeff,curTime, timeFinal))*sin(calc_gamma(coeff,curTime, timeFinal)));
	
}

template <class T> T calcRate_vtheta(const elements<T> & y, coefficients<T> & coeff, const T & accel, const T & curTime, const T & timeFinal)
{
	return -y.vr*y.vtheta / y.r + accel*cos(calc_tau(coeff,curTime, timeFinal))*cos(calc_gamma(coeff,curTime, timeFinal));

}

template <class T> T calcRate_vz(const elements<T> & y, coefficients<T> & coeff, const T & accel, const T & curTime, const T & timeFinal)
{
	return (-constG * massSun * y.z / (pow(pow(y.r, 2) + pow(y.z, 2), 1.5))) + accel*sin(calc_tau(coeff,curTime, timeFinal));
	
}

//////////////////////////////////////////////////////////////////////////
// Orbital motion equations - only acceleration changes
//////////////////////////////////////////////////////////////////////////



template <class T> T calcRate_vr_earth(const elements<T> & y, const T & curTime, const T & timeFinal)
{
	return (-constG * massSun * y.r / (pow(pow(y.r, 2) + pow(y.z, 2), 1.5))) + (pow(y.vtheta,2) / y.r); 
	
}

template <class T> T calcRate_vtheta_earth(const elements<T> & y, const T & curTime, const T & timeFinal)
{
	return -y.vr*y.vtheta / y.r;
}

template <class T> T calcRate_vz_earth(const elements<T> & y, const T & curTime, const T & timeFinal)
{
	return (-constG * massSun * y.z / pow(pow(y.r, 2) + pow(y.z, 2), 1.5));
	
}