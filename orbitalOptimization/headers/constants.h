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

// starting location in the optimization array
#define GAMMA_OFFSET 0 
#define TAU_OFFSET 9
#define ALPHA_OFFSET 12
#define BETA_OFFSET 13 
#define TRIPTIME_OFFSET 14
#define COAST_OFFSET 15
#define COASTTHRESHOLD_OFFSET 18

// Spacecraft constants
#define wetMass 3000 // fuel + dry mass (kg)
#define dryMass 2900 // (kg)

#endif