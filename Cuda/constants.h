#ifndef constants_h
#define constants_h
#include <math.h>

#define _USE_MATH_DEFINES // for use of M_PI

// Planetary properties and constants
#define earthRadius 1.49598261e11/AU // radial distance of Earth from Sun (au)
#define earthMass 5.97237e24 // mass of the earth (kg)
#define lunarMass 7.34767e22 // mass of the earth's moon (kg)
#define ESOI earthRadius*pow(((earthMass+lunarMass)/massSun),0.4) // escape sphere of influence (au)
#define AU 1.49597870691e11 // used to convert meters to astronomical units (m) 
#define constG 1.99349603314131e-44 // gravitational constant- used to calculate the gravitational force (AU^3/(s^2 * kg)) 
#define massSun 1.988500e30 // mass of the sun (kg)
#define orbitalPeriod 6.653820100923719e+07 // orbital period time of the asteroid (s)
#define orbitalInclination 0.0594906 //orbital inclination of the asteroid (rad)


// Starting location and sizes in the optimization array for navigation to access specific values
#define OPTIM_VARS 19 // Number of variables in the optimization
#define GAMMA_ARRAY_SIZE 7 // Length of the array of coefficients for gamma
#define TAU_ARRAY_SIZE 3 // Length of the array of coefficients for tau
#define COAST_ARRAY_SIZE 5 // Length of the array of coefficients for coasting

#define GAMMA_OFFSET 0 // x[0-6] fourth order fourier for in-plane angle
#define TAU_OFFSET 7 // x[7-9] first order fourier for out-of-plane angle
#define ALPHA_OFFSET 10 // x[10] position escape earth angle
#define BETA_OFFSET 11 // x[11] velocity escape in-plane earth angle
#define ZETA_OFFSET 12 // x[12] velocity escape out-of-plane earth angle
#define TRIPTIME_OFFSET 13 // x[13] total duration of the trip
#define COAST_OFFSET 14 // x[14-16] second order fourier for coasting determination

#define ASTEROID_OFFSET 0 // y0A
#define EARTH_OFFSET 6 // y0E
#define THRESHOLD_OFFSET 12 // coast threshold

#endif