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
#define gamma_offset 0 
#define tau_offset 7 
#define alpha_offset 10
#define beta_offset 11 
#define tripTime_offset 12 

// Spacecraft constants
#define wetMass 3000 // fuel + dry mass (kg)
#define NEXTmDot 5.76*pow(kiloConversion,2)// fuel differential (kg/s)
#define NEXTeff 0.671 // NEXT efficiency (unitless)
#define NEXTP0 7.240/kiloConversion // initial power (W)

#endif