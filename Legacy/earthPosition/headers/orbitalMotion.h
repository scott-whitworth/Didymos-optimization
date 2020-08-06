#ifndef orbitalMotion_h
#define orbitalMotion_h
#include "elements.h"

// Input for earthInitial:
// trip time - total time from launch to impact, sets the initial earth position
// calculates the earth's initial conditions for launch date based on impact date (oct. 5, 2022) minus the trip time
elements<double> earthInitial(double timeInitial, double tripTime,const elements<double> & earth);



#include "orbitalMotion.cpp"
#endif