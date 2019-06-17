#ifndef constants_h
#define constants_h

#define earthRadius 1.49598261e11/AU // (au)
#define earthMass 5.9742e24 // (kg)
#define ESOI earthRadius*pow((earthMass/massSun),0.4) // escape sphere of influence (au)
#define C3Energy 4.676e6 // (m^2/s^2)
#define vEscape sqrt(2*C3Energy)/AU // (au/s)
#define AU 1.49597870691e11// units: m; used to convert meters to astronomical units
#define constG 1.99349603314131e-44 //units: AU^3/(s^2 * kg); gravitational constant- used to calculate the gravitational force
#define massSun 1.98892e30//kg
#define Torbital 6.653820100923719e+07 // orbital period time of the asteroid (s)

#endif