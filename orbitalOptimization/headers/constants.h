#ifndef constants_h
#define constants_h
#include <math.h>

#define earthRadius 1.49598261e11/AU // radius of the earth (au)
#define earthMass 5.9742e24 // mass of the earth (kg)
#define ESOI earthRadius*pow((earthMass/massSun),0.4) // escape sphere of influence (au)
#define C3Energy 4.676e6 // energy of spacecraft at earth escape (m^2/s^2)
#define vEscape sqrt(C3Energy)/AU // magnitude of velocity at earth escape  (au/s)
#define AU 1.49597870691e11// used to convert meters to astronomical units (m) 
#define constG 1.99349603314131e-44 //gravitational constant- used to calculate the gravitational force (AU^3/(s^2 * kg)) 
#define massSun 1.98892e30// mass of the sun (kg)
#define Torbital 6.653820100923719e+07 // orbital period time of the asteroid (s)
#define kiloConversion .001 //  used to convert to kilograms and Watts (meters/kilometer)


#define MAX_NUMSTEPS 1e9
#define RK_TOL 1e-12
#define F_MIN 1e-20

//final conditions of the asteroid and Earth
#define R_FIN_AST 1.02696822710421
#define THETA_FIN_AST 0.238839574416454
#define Z_FIN_AST -0.0526614832914496
#define VR_FIN_AST -2.05295246185041e-08
#define VTHETA_FIN_AST 2.29132593453064e-07
#define VZ_FIN_AST 8.00663905822009e-09

#define R_FIN_EARTH 1.00021392223428
#define THETA_FIN_EARTH 0.199470650149394
#define Z_FIN_EARTH -1.54878511585620e-05
#define VR_FIN_EARTH -3.32034068725821e-09
#define VTHETA_FIN_EARTH 1.99029138292504e-07
#define VZ_FIN_EARTH -9.71518257891386e-12

//Official DART mission data
#define V_IMPACT 4.399378072e-08 //impact velocity in AU/s
//#define V_IMPACT 6.0e-08 //impact velocity in AU/s

// starting location in the optimization array
#define OPTIM_VARS 21//Number of variables in the optimization
#define GAMMA_OFFSET 0 // x[0-8] fourth order fourier for in-plane angle
#define TAU_OFFSET 9 // x[9-11] first order fourier for out-of-plane angle
#define ALPHA_OFFSET 12 // x[12] position escape earth angle
#define BETA_OFFSET 13 // x[13] velocity escape earth angle
#define TRIPTIME_OFFSET 14 // x[14] total duration of the trip
#define COAST_OFFSET 15 // x[15-19] second order fourier for coasting determination
#define THRESHOLD_OFFSET 20 // x[20] coasting threshold
//#define WETMASS_OFFSET 21 // x[21] used to calculate the amount of fuel avaliable to the spacecraft

// Spacecraft constants
#define dryMass 2700 // mass of spacecraft excluding fuel (kg)
#define WET_MASS 3000

#endif