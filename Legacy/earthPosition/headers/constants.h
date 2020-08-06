#ifndef constants_h
#define constants_h

// Planetary Constants
#define constG 1.99349603314131e-44 //gravitational constant- used to calculate the gravitational force (AU^3/(s^2 * kg)) 
#define AU 1.49597870691e11// used to convert meters to astronomical units (m) 
#define earthRadius 1.49598261e11/AU // radius of the earth (au)
#define earthMass 5.9742e24 // mass of the earth (kg)
#define massSun 1.98892e30// mass of the sun (kg)

// Constants used in functions
#define MAX_NUMSTEPS 1e9 // The highest precision the runge kutta method is going to use for the first step
#define RK_TOL 1e-12 // The relative/absolute (not sure which one it is) tolerance for the runge kutta algorithm

// Final conditions of the earth on impact date (oct 5, 2022)
#define R_FIN_EARTH 1.00021392223428 // radial position (au)
#define THETA_FIN_EARTH 0.199470650149394 // angular postion (rad)
#define Z_FIN_EARTH -1.54878511585620e-05 // axial position (au)
#define VR_FIN_EARTH -3.32034068725821e-09 // radial velocity (au/s)
#define VTHETA_FIN_EARTH 1.99029138292504e-07 // angular velocity (au/s)
#define VZ_FIN_EARTH -9.71518257891386e-12 // axial velocity (au/s)

#endif